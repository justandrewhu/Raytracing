#pragma once

#include "scene.hpp"
#include "camera.hpp"
#include "sphere.hpp"
#include "triangle.hpp"
#include "mesh.hpp"
#include "csg.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <stdexcept>

namespace rt {

// Simple JSON-like parser for scene files
// Format is a subset of JSON, supporting objects, arrays, strings, numbers
class SimpleJSONParser {
public:
    enum class TokenType { LBRACE, RBRACE, LBRACKET, RBRACKET, COLON, COMMA, STRING, NUMBER, TRUE, FALSE, NULL_VAL, END };

    struct Token {
        TokenType type;
        std::string value;
    };

    struct Value {
        enum class Type { Null, Bool, Number, String, Array, Object };
        Type type = Type::Null;
        bool bool_val = false;
        double num_val = 0;
        std::string str_val;
        std::vector<Value> array_val;
        std::map<std::string, Value> object_val;

        bool is_null() const { return type == Type::Null; }
        bool is_bool() const { return type == Type::Bool; }
        bool is_number() const { return type == Type::Number; }
        bool is_string() const { return type == Type::String; }
        bool is_array() const { return type == Type::Array; }
        bool is_object() const { return type == Type::Object; }

        double as_number(double default_val = 0) const {
            return is_number() ? num_val : default_val;
        }

        std::string as_string(const std::string& default_val = "") const {
            return is_string() ? str_val : default_val;
        }

        bool as_bool(bool default_val = false) const {
            return is_bool() ? bool_val : default_val;
        }

        const Value& operator[](const std::string& key) const {
            static Value null_val;
            if (!is_object()) return null_val;
            auto it = object_val.find(key);
            return it != object_val.end() ? it->second : null_val;
        }

        const Value& operator[](size_t index) const {
            static Value null_val;
            if (!is_array() || index >= array_val.size()) return null_val;
            return array_val[index];
        }

        size_t size() const {
            if (is_array()) return array_val.size();
            if (is_object()) return object_val.size();
            return 0;
        }

        bool has(const std::string& key) const {
            return is_object() && object_val.find(key) != object_val.end();
        }
    };

private:
    std::string input;
    size_t pos = 0;

    char peek() { return pos < input.size() ? input[pos] : '\0'; }
    char get() { return pos < input.size() ? input[pos++] : '\0'; }
    void skip_whitespace() { while (pos < input.size() && std::isspace(input[pos])) ++pos; }

    Token next_token() {
        skip_whitespace();
        if (pos >= input.size()) return {TokenType::END, ""};

        char c = peek();

        if (c == '{') { get(); return {TokenType::LBRACE, "{"}; }
        if (c == '}') { get(); return {TokenType::RBRACE, "}"}; }
        if (c == '[') { get(); return {TokenType::LBRACKET, "["}; }
        if (c == ']') { get(); return {TokenType::RBRACKET, "]"}; }
        if (c == ':') { get(); return {TokenType::COLON, ":"}; }
        if (c == ',') { get(); return {TokenType::COMMA, ","}; }

        if (c == '"') {
            get();
            std::string str;
            while (peek() != '"' && peek() != '\0') {
                if (peek() == '\\') {
                    get();
                    char esc = get();
                    switch (esc) {
                        case 'n': str += '\n'; break;
                        case 't': str += '\t'; break;
                        case 'r': str += '\r'; break;
                        case '"': str += '"'; break;
                        case '\\': str += '\\'; break;
                        default: str += esc;
                    }
                } else {
                    str += get();
                }
            }
            get(); // closing quote
            return {TokenType::STRING, str};
        }

        if (c == '-' || std::isdigit(c)) {
            std::string num;
            if (c == '-') num += get();
            while (std::isdigit(peek())) num += get();
            if (peek() == '.') {
                num += get();
                while (std::isdigit(peek())) num += get();
            }
            if (peek() == 'e' || peek() == 'E') {
                num += get();
                if (peek() == '+' || peek() == '-') num += get();
                while (std::isdigit(peek())) num += get();
            }
            return {TokenType::NUMBER, num};
        }

        if (std::isalpha(c)) {
            std::string word;
            while (std::isalpha(peek())) word += get();
            if (word == "true") return {TokenType::TRUE, "true"};
            if (word == "false") return {TokenType::FALSE, "false"};
            if (word == "null") return {TokenType::NULL_VAL, "null"};
            throw std::runtime_error("Unknown keyword: " + word);
        }

        throw std::runtime_error(std::string("Unexpected character: ") + c);
    }

    Token current;
    Token lookahead;

    void advance() {
        current = lookahead;
        lookahead = next_token();
    }

    void expect(TokenType type) {
        if (current.type != type) {
            throw std::runtime_error("Unexpected token");
        }
        advance();
    }

    Value parse_value() {
        Value v;
        switch (current.type) {
            case TokenType::NULL_VAL:
                v.type = Value::Type::Null;
                advance();
                break;
            case TokenType::TRUE:
                v.type = Value::Type::Bool;
                v.bool_val = true;
                advance();
                break;
            case TokenType::FALSE:
                v.type = Value::Type::Bool;
                v.bool_val = false;
                advance();
                break;
            case TokenType::NUMBER:
                v.type = Value::Type::Number;
                v.num_val = std::stod(current.value);
                advance();
                break;
            case TokenType::STRING:
                v.type = Value::Type::String;
                v.str_val = current.value;
                advance();
                break;
            case TokenType::LBRACKET:
                v.type = Value::Type::Array;
                advance();
                if (current.type != TokenType::RBRACKET) {
                    v.array_val.push_back(parse_value());
                    while (current.type == TokenType::COMMA) {
                        advance();
                        v.array_val.push_back(parse_value());
                    }
                }
                expect(TokenType::RBRACKET);
                break;
            case TokenType::LBRACE:
                v.type = Value::Type::Object;
                advance();
                if (current.type != TokenType::RBRACE) {
                    std::string key = current.value;
                    expect(TokenType::STRING);
                    expect(TokenType::COLON);
                    v.object_val[key] = parse_value();
                    while (current.type == TokenType::COMMA) {
                        advance();
                        key = current.value;
                        expect(TokenType::STRING);
                        expect(TokenType::COLON);
                        v.object_val[key] = parse_value();
                    }
                }
                expect(TokenType::RBRACE);
                break;
            default:
                throw std::runtime_error("Unexpected token in value");
        }
        return v;
    }

public:
    Value parse(const std::string& json) {
        input = json;
        pos = 0;
        lookahead = next_token();
        advance();
        return parse_value();
    }
};


class SceneLoader {
public:
    using JSON = SimpleJSONParser::Value;

    struct SceneData {
        Camera camera;
        Scene scene;
    };

    static Vec3 parse_vec3(const JSON& j, const Vec3& default_val = Vec3(0)) {
        if (j.is_array() && j.size() >= 3) {
            return Vec3(j[0].as_number(), j[1].as_number(), j[2].as_number());
        }
        if (j.is_number()) {
            return Vec3(j.as_number());
        }
        return default_val;
    }

    static Color parse_color(const JSON& j, const Color& default_val = Color(0)) {
        return parse_vec3(j, default_val);
    }

    static MaterialPtr parse_material(const JSON& j) {
        auto mat = std::make_shared<Material>();

        mat->k_d = parse_color(j["diffuse"], Color(0.5));
        mat->k_s = parse_color(j["specular"], Color(0));
        mat->p = j["exponent"].as_number(20.0);
        mat->k_m = parse_color(j["mirror"], Color(0));
        mat->k_a = j.has("ambient") ? parse_color(j["ambient"]) : mat->k_d;
        mat->ior = j["ior"].as_number(1.0);
        mat->k_t = parse_color(j["transmission"], Color(0));
        mat->k_e = parse_color(j["emission"], Color(0));

        return mat;
    }

    static std::string get_directory(const std::string& filepath) {
        size_t pos = filepath.find_last_of("/\\");
        return (pos == std::string::npos) ? "." : filepath.substr(0, pos);
    }

    static std::string resolve_path(const std::string& path, const std::string& base_dir) {
        if (!path.empty() && (path[0] == '/' || (path.size() > 1 && path[1] == ':'))) {
            return path;
        }
        // Relative path - resolve from base directory
        return base_dir + "/" + path;
    }

    static ShapePtr parse_shape(const JSON& j, std::map<std::string, MaterialPtr>& materials,
                                const std::string& base_dir = ".") {
        std::string type = j["type"].as_string();
        MaterialPtr mat = nullptr;

        if (j.has("material")) {
            if (j["material"].is_string()) {
                mat = materials[j["material"].as_string()];
            } else {
                mat = parse_material(j["material"]);
            }
        }

        if (!mat) {
            mat = std::make_shared<Material>();
        }

        if (type == "sphere") {
            Vec3 center = parse_vec3(j["center"]);
            double radius = j["radius"].as_number(1.0);
            return std::make_shared<Sphere>(center, radius, mat);
        }

        if (type == "triangle") {
            Vec3 v0 = parse_vec3(j["v0"]);
            Vec3 v1 = parse_vec3(j["v1"]);
            Vec3 v2 = parse_vec3(j["v2"]);
            return std::make_shared<Triangle>(v0, v1, v2, mat);
        }

        if (type == "mesh") {
            std::string filename = j["file"].as_string();
            std::string resolved_path = resolve_path(filename, base_dir);
            auto mesh = Mesh::load_obj(resolved_path, mat);

            // Parse transformations
            Vec3 scale_vec(1, 1, 1);
            Vec3 rotation_vec(0, 0, 0);
            Vec3 translation_vec(0, 0, 0);

            if (j.has("scale")) {
                if (j["scale"].is_number()) {
                    double s = j["scale"].as_number(1.0);
                    scale_vec = Vec3(s, s, s);
                } else {
                    scale_vec = parse_vec3(j["scale"], Vec3(1, 1, 1));
                }
            }

            if (j.has("rotation")) {
                rotation_vec = parse_vec3(j["rotation"], Vec3(0, 0, 0));
            }

            if (j.has("translation")) {
                translation_vec = parse_vec3(j["translation"], Vec3(0, 0, 0));
            }

            // Apply transformations: scale -> rotate -> translate
            Mat3 rot_matrix = Mat3::rotate(rotation_vec.x, rotation_vec.y, rotation_vec.z);
            mesh->transform(rot_matrix, scale_vec, translation_vec);

            return mesh;
        }

        if (type == "box") {
            Vec3 min_corner = parse_vec3(j["min"]);
            Vec3 max_corner = parse_vec3(j["max"]);
            return make_box(min_corner, max_corner, mat);
        }

        if (type == "csg_union" || type == "csg_intersection" || type == "csg_difference") {
            ShapePtr left = parse_shape(j["left"], materials, base_dir);
            ShapePtr right = parse_shape(j["right"], materials, base_dir);

            if (type == "csg_union") {
                return csg_union(left, right, mat);
            } else if (type == "csg_intersection") {
                return csg_intersection(left, right, mat);
            } else {
                return csg_difference(left, right, mat);
            }
        }

        throw std::runtime_error("Unknown shape type: " + type);
    }

    static LightPtr parse_light(const JSON& j) {
        std::string type = j["type"].as_string();

        if (type == "point") {
            Vec3 position = parse_vec3(j["position"]);
            Color intensity = parse_color(j["intensity"], Color(100));
            return std::make_shared<PointLight>(position, intensity);
        }

        if (type == "ambient") {
            Color intensity = parse_color(j["intensity"], Color(0.1));
            return std::make_shared<AmbientLight>(intensity);
        }

        throw std::runtime_error("Unknown light type: " + type);
    }

    static Camera parse_camera(const JSON& j) {
        Vec3 eye = parse_vec3(j["eye"], Vec3(0, 0, 5));
        Vec3 target = parse_vec3(j["target"], Vec3(0, 0, 0));
        Vec3 up = parse_vec3(j["up"], Vec3(0, 1, 0));
        double vfov = j["vfov"].as_number(45.0);
        double aspect = j["aspect"].as_number(16.0 / 9.0);

        return Camera(eye, target, up, vfov, aspect);
    }

    static SceneData load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open scene file: " + filename);
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string content = buffer.str();

        SimpleJSONParser parser;
        JSON root = parser.parse(content);

        std::string base_dir = get_directory(filename);

        SceneData data;

        if (root.has("background")) {
            data.scene.background_color = parse_color(root["background"], Color(0.2, 0.3, 0.5));
        }

        std::map<std::string, MaterialPtr> materials;
        if (root.has("materials") && root["materials"].is_object()) {
            for (const auto& [name, mat_json] : root["materials"].object_val) {
                materials[name] = parse_material(mat_json);
            }
        }

        if (root.has("shapes") && root["shapes"].is_array()) {
            for (size_t i = 0; i < root["shapes"].size(); ++i) {
                data.scene.add_shape(parse_shape(root["shapes"][i], materials, base_dir));
            }
        }

        if (root.has("lights") && root["lights"].is_array()) {
            for (size_t i = 0; i < root["lights"].size(); ++i) {
                data.scene.add_light(parse_light(root["lights"][i]));
            }
        }

        if (root.has("camera")) {
            data.camera = parse_camera(root["camera"]);
        }

        // speeedy
        data.scene.build_bvh();

        return data;
    }
};

} // namespace rt
