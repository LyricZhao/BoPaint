# include "graph.h"

# include <nlohmann/json.hpp>
# include <opencv2/opencv.hpp>

# include <cassert>
# include <cmath>
# include <fstream>
# include <iostream>

using json = nlohmann:: json;

const double sampling_weigh[3][3] = {
    1 / 16., 2 / 16., 1 / 16.,
    2 / 16., 4 / 16., 2 / 16.,
    1 / 16., 2 / 16., 1 / 16.
};

const double eps = 1e-6;

inline void debug() {
    std:: cout << "DEBUG" << std:: endl;
    return;
}

Graph:: Graph(std:: string input_path) {
    img = nullptr;
    rev_xy = false;
    waaa = false, elaa = 1;
    height = width = 0;
    readFromFile(input_path);
    return;
}

Graph:: ~Graph() {
    if(img != nullptr) {
        delete img;
    }
    return;
}

void Graph:: readFromFile(std:: string path) {
    std:: ifstream input(path);
    json config; input >> config;
    height = config["height"], width = config["width"];
    waaa = config["waaa"], elaa = config["elaa"];
    for(auto line: config["lines"]) {
        Point a = Point(line[0], line[1]);
        Point b = Point(line[2], line[3]);
        lines.push_back(Line(a, b));
    }
    for(auto arc: config["arcs"]) {
        Point a = Point(arc[0], arc[1]);
        Point b = Point(arc[2], arc[3]);
        arcs.push_back(Arc(a, b, arc[4]));
    }
    return;
}

void Graph:: show() {
    cv:: imshow("output", *img);
    cv:: waitKey(0);
    cv:: destroyAllWindows();
    return;
}

void Graph:: saveToFile(std:: string path) {
    cv:: imwrite(path, *img);
    return;
}

bool Graph:: inRange(Point a) {
    return 0 <= a.x && a.x <= height - 1 && 0 <= a.y && a.y <= width - 1;
}

void Graph:: drawPixel(Point a, Color color) {
    if(rev_xy) std:: swap(a.x, a.y);
    img -> at<cv:: Vec3b>(a.x, a.y)[0] = (ColorBit) ((color & 0x00ff0000ul) >> 16); /* B */
    img -> at<cv:: Vec3b>(a.x, a.y)[1] = (ColorBit) ((color & 0x0000ff00ul) >>  8); /* G */
    img -> at<cv:: Vec3b>(a.x, a.y)[2] = (ColorBit) ((color & 0x000000fful) >>  0); /* R */
    return;
}

void Graph:: drawPixel(Point a, ColorBit b, ColorBit g, ColorBit r) {
    if(rev_xy) std:: swap(a.x, a.y);
    img -> at<cv:: Vec3b>(a.x, a.y)[0] = b;
    img -> at<cv:: Vec3b>(a.x, a.y)[1] = g;
    img -> at<cv:: Vec3b>(a.x, a.y)[2] = r;
    return;
}

void Graph:: weighedSampling(double lk, double lb, Color color, Point a) {
    if(!inRange(a)) return;
    double gray = 0;
    double ulb = lb + sqrt(1 + lk * lk) / 2.;
    double dlb = lb - sqrt(1 + lk * lk) / 2.;
    for(int i = 0; i < 3; ++ i) {
        for(int j = 0; j < 3; ++ j) {
            double lux = a.x + i * 1 / 3., luy = a.y + (j + 1) * 1 / 3.;
            double ldx = a.x + i * 1 / 3., ldy = a.y + j * 1 / 3.;
            double rdx = a.x + (i + 1) * 1 / 3., rdy = a.y + j * 1 / 3.;
            double rux = a.x + (i + 1) * 1 / 3., ruy = a.y + (j + 1) * 1 / 3.;
            if((lk >= 0 && (lk * rdx + ulb > rdy + eps && lk * lux + dlb < luy - eps)) ||
                (lk < 0 && (lk * ldx + ulb > ldy + eps && lk * rux + dlb < ruy - eps)) ) {
                gray += sampling_weigh[i][j];
            }
        }
    }
    if(gray > eps) {
        drawPixel(a, (ColorBit) ((255 - (255 - ((color & 0x00ff0000ul) >> 16)) * gray)),
                     (ColorBit) ((255 - (255 - ((color & 0x0000ff00ul) >>  8)) * gray)),
                     (ColorBit) ((255 - (255 - ((color & 0x000000fful) >>  0)) * gray)));
    }
    return;
}

void Graph:: lineAA(Point a, double lk, double lb, Color color) {
    Point shift(0, 1);
    weighedSampling(lk, lb, color, a);
    weighedSampling(lk, lb, color, a - shift);
    weighedSampling(lk, lb, color, a + shift);
    return;
}

void Graph:: bresenhamLine(Line line) {
    int dx, dy, e;
    double lk, lb;
    Point a = line.first, b = line.second;
    if(fabs(a.x - b.x) < fabs(a.y - b.y)) {
        rev_xy = true;
        std:: swap(a.x, a.y);
        std:: swap(b.x, b.y);
    }
    if(a.x > b.x) std:: swap(a, b);
    dx = b.x - a.x, dy = b.y - a.y, e = -dx;
    lk = 1.0 * dy / dx, lb = a.y - lk * a.x;
    for(int i = 0; i <= dx; ++ i) {
        if(waaa) lineAA(a, lk, lb, BLACK);
        else drawPixel(a, BLACK);
        a.x += 1, e += 2 * dy;
        if(e >= 0) {
            a.y += 1, e -= 2 * dx;
        }
    }
    if(rev_xy) rev_xy = false;
    return;
}

void Graph:: draw() {
    img = new cv:: Mat(height, width, CV_8UC3, cv:: Scalar(255, 255, 255));
    for(auto line: lines) {
        bresenhamLine(line);
    }
    for(auto arc: arcs) {

    }
    return;
}