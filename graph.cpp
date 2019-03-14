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
const double pi = 3.14159265;

inline void debug() {
    std:: cout << "DEBUG" << std:: endl;
    return;
}

inline bool inRange(double deg, Range range) {
    return range.first - eps <= deg && deg <= range.second + eps;
}

inline double sqr(double x) {
    return x * x;
}

inline bool inCircle(double r, double vx, double vy) {
    return sqr(vx) + sqr(vy) <= sqr(r);
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
        Point center = Point(arc[0], arc[1]);
        Range range = Range(arc[3], arc[4]);
        arcs.push_back(Arc(center, arc[2], range));
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
    if(!inRange(a)) return;
    if(rev_xy) std:: swap(a.x, a.y);
    img -> at<cv:: Vec3b>(a.x, a.y)[0] = (ColorBit) ((color & 0x00ff0000ul) >> 16); /* B */
    img -> at<cv:: Vec3b>(a.x, a.y)[1] = (ColorBit) ((color & 0x0000ff00ul) >>  8); /* G */
    img -> at<cv:: Vec3b>(a.x, a.y)[2] = (ColorBit) ((color & 0x000000fful) >>  0); /* R */
    return;
}

void Graph:: drawPixel(Point a, ColorBit b, ColorBit g, ColorBit r) {
    if(!inRange(a)) return;
    if(rev_xy) std:: swap(a.x, a.y);
    img -> at<cv:: Vec3b>(a.x, a.y)[0] = b;
    img -> at<cv:: Vec3b>(a.x, a.y)[1] = g;
    img -> at<cv:: Vec3b>(a.x, a.y)[2] = r;
    return;
}

void Graph:: drawPixel(Point a, Color color, double gray) {
    if(gray < eps || !inRange(a)) return;
    drawPixel(a, (ColorBit) ((255 - (255 - ((color & 0x00ff0000ul) >> 16)) * gray)),
                              (ColorBit) ((255 - (255 - ((color & 0x0000ff00ul) >>  8)) * gray)),
                              (ColorBit) ((255 - (255 - ((color & 0x000000fful) >>  0)) * gray)));
    return;
}

void Graph:: weighedSamplingLine(double lk, double lb, Color color, Point a) {
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
    drawPixel(a, color, gray);
    return;
}

double Graph:: weighedSamplingCircle(int r, Point delta) {
    double gray = 0;
    double lr = r - 0.3, rr = r + 0.3;
    for(int i = 0; i < 3; ++ i) {
        for(int j = 0; j < 3; ++ j) {
            bool flag = false;
            for(int k = 0; k < 4 && !flag; ++ k) {
                double x = delta.x + (i + (k & 1)) * 1 / 3.;
                double y = delta.y + (i + ((k >> 1) & 1)) * 1 / 3.;
                flag |= !inCircle(lr, x, y) && !inCircle(rr, x, y);
            }
            if(flag) {
                gray += sampling_weigh[i][j];
            }
        }
    }
    return gray;
}

void Graph:: lineAA(Point a, double lk, double lb, Color color) {
    Point shift(0, 1);
    weighedSamplingLine(lk, lb, color, a);
    weighedSamplingLine(lk, lb, color, a - shift);
    weighedSamplingLine(lk, lb, color, a + shift);
    return;
}

void Graph:: bresenhamLine(Line line) {
    int dx, dy, e, sign_y;
    double lk, lb;
    Point a = line.first, b = line.second;
    assert(!(a == b));
    if(fabs(a.x - b.x) < fabs(a.y - b.y)) {
        rev_xy = true;
        std:: swap(a.x, a.y);
        std:: swap(b.x, b.y);
    }
    if(a.x > b.x) std:: swap(a, b);
    dx = b.x - a.x, dy = b.y - a.y;
    sign_y = (dy < 0) ? -1 : 1;
    e = -sign_y * dx;
    lk = 1.0 * dy / dx, lb = a.y - lk * a.x;
    for(int i = 0; i <= dx; ++ i) {
        if(waaa) lineAA(a, lk, lb, BLACK);
        else drawPixel(a, BLACK);
        a.x += 1, e += 2 * dy;
        if(e * sign_y >= 0) {
            a.y += sign_y, e += -2 * sign_y * dx;
        }
    }
    if(rev_xy) rev_xy = false;
    return;
}

void Graph:: circle8Point(Point delta, Point center, int r, Range range, Color color) {
    const int sign[2] = {-1, 1};
    Point shift = Point(0, 1), eta[3];
    Point aa_delta[3] = {delta, delta - shift, delta + shift};
    double aa_gray[3];
    for(int i = 0; i < (waaa ? 3 : 0); ++ i) {
        aa_gray[i] = weighedSamplingCircle(r, aa_delta[i]);
    }
    for(int i = 0; i < 2; ++ i) {
        for(int j = 0; j < 2; ++ j) {
            for(int k = 0; k < 2; ++ k) {
                for(int l = 0; l < 3; ++ l) {
                    eta[l] = aa_delta[l];
                    if(i) std:: swap(eta[l].x, eta[l].y);
                    eta[l].x *= sign[j], eta[l].y *= sign[k];
                }
                if(:: inRange(atan2(eta[0].y, eta[0].x) * 180 / pi, range)) {
                    if(waaa) {
                        drawPixel(center + eta[0], color, aa_gray[0]);
                        drawPixel(center + eta[1], color, aa_gray[1]);
                        drawPixel(center + eta[2], color, aa_gray[2]);
                    } else {
                        drawPixel(center + eta[0], color);
                    }
                }
            }
        }
    }
    return;
}

void Graph:: drawCircle(Point center, int r, Range range, Color color) {
    Point delta = Point(0, r);
    double d = 1.25 - r;
    circle8Point(delta, center, r, range, color);
    while(delta.x <= delta.y){
        if(d < 0) {
            d += 2. * delta.x + 3.;
        }
        else {
            d += 2 * (delta.x - delta.y) + 5;
            delta.y -= 1;
        }
        delta.x += 1;
        circle8Point(delta, center, r, range, color);
    }
    return;
}

void Graph:: draw() {
    img = new cv:: Mat(height, width, CV_8UC3, cv:: Scalar(255, 255, 255));
    for(auto line: lines) {
        bresenhamLine(line);
    }
    for(auto arc: arcs) {
        drawCircle(std:: get<0>(arc), std:: get<1>(arc), std:: get<2>(arc), BLACK);
    }
    return;
}