// Author: Jiahui Huang from CSCG.

#include <iostream>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;

int main(int argc, char** argv) {
	
	// Parameters.
	const int width = 300;
	const int height = 200;

	// Build a canvas. Each element in the canvas is 3 unsigned char (8UC3) for in the order of BGR.
	cv::Mat img(height, width, CV_8UC3, cv::Scalar(0, 0, 0));

	// Draw an antialiased blue line. Don't use this function in your homework.
	cv::line(img, cv::Point(10, 10), cv::Point(100, 100), cv::Scalar(255, 0, 0), 1, LINE_AA);

	// Draw a horizontal line pixel by pixel.
	// Note that this can be done more efficiently by using techniques from:
	// https://docs.opencv.org/2.4/doc/tutorials/core/how_to_scan_images/how_to_scan_images.html
	for (int c = 0; c < width; ++c) {
		img.at<cv::Vec3b>(100, c)[1] = 255;
	}

	// Save the image.
	cv::imwrite("output.png", img);

	// Show the image, wait for user keystroke and quit.
	cv::imshow("output", img);
	cv::waitKey(0);
	cv::destroyAllWindows();
	return 0;
}
