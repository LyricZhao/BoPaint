# include <nlohmann/json.hpp>
# include <opencv2/opencv.hpp>

# include <algorithm>
# include <string>
# include <tuple>
# include <vector>

# define BLACK 0x00000000ul

using json = nlohmann:: json;

struct Point {
    int x, y;
    Point(int _x, int _y): x(_x), y(_y) {}
    friend Point operator + (const Point &a, const Point &b) { return Point(a.x + b.x, a.y + b.y); }
    friend Point operator - (const Point &a, const Point &b) { return Point(a.x - b.x, a.y - b.y); }
    friend bool operator < (const Point &a, const Point &b) {
        return a.x == b.x ? a.y < b.y : a.x < b.x;
    }
    friend bool operator > (const Point &a, const Point &b) {
        return a.x == b.x ? a.y > b.y : a.x > b.x;
    }
    friend bool operator == (const Point &a, const Point &b) {
        return a.x == b.x && a.y == b.y;
    }
    void print() {
        std:: cout << "(" << x << ", " << y << ")";
        return;
    }
};

using Arc = std:: tuple<Point, Point, int>;
using Line = std:: pair<Point, Point>;

typedef unsigned int Color;
typedef unsigned char ColorBit;

class Graph {
public:
    Graph(std:: string input_path);
    ~Graph();
    void draw();
    void show();
    void readFromFile(std:: string path);
    void saveToFile(std:: string path);

private:
    bool rev_xy;
    int height, width;
    cv:: Mat *img;
    bool waaa; int elaa;
    std:: vector<Arc> arcs;
    std:: vector<Line> lines;

    bool inRange(Point a);
    void lineAA(Point a, double lk, double lb, Color color);
    void weighedSampling(double lk, double lb, Color color, Point a);
    void drawPixel(Point a, Color color);
    void drawPixel(Point a, ColorBit r, ColorBit g, ColorBit b);
    void bresenhamLine(Line line);
};