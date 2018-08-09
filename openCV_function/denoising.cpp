/***
 * Must be compiled with:
 * g++ denoising.cpp `pkg-config --cflags --libs opencv`
 * 
 ***/

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/photo/photo.hpp"
#include "opencv2/opencv.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>

using namespace cv;
using namespace std;

cv::Mat src, src_gray, denoised, detectedEdges, bw;
char targetName [20] = "denoised";
char srcName[20] = "original";

int main(int argc, char** argv) 
{
  src = cv::imread(argv[1]);
	//src = cv::imread("102.tif");
  if( !src.data )
  {
    std::cout << "no img could be read" << std::endl;
  }

  cv::cvtColor(src, src_gray, CV_BGR2GRAY);
 
  /*
   * opencv built-in function for denoising
   *
   * src: 8-bit image
   * target: must be similar to src
   * filter strength: float
   * templateWindowSize: int
   * searchWindowSize: int
   *
   */
  cv::fastNlMeansDenoising(src_gray, denoised,  5, 23, 31);
  std::stringstream ss;
  ss << "denoised_b_"<< argv[1] ;
  cv::imwrite(ss.str(), denoised);

  bw.create(denoised.size(), denoised.type());
  //cv::adaptiveThreshold(denoised, bw, 255, ADAPTIVE_THRESH_MEAN_C,
  //    THRESH_BINARY, 31, 0);
  cv::threshold(denoised, bw, 0, 255, CV_THRESH_OTSU);
  std::stringstream ss2;
  ss2 << "binarized_b_"<< argv[1];
  cv::imwrite(ss2.str(), bw);
  /*
   * Use the denoised image for threshold and edge detection
   */
  /*
  //for Canny filter
  Mat dst;
  int kernel = 3;
  int edgeThresh = 1;
  int lowThresh = 0;
  int const max_lowThresh=100;
  int ratio = 3;
  //for contours
  int thresh = 100;
  int max_thresh = 255;
  RNG rng(12345);

  denoised.create(src_gray.size(), src_gray.type()); 
  detectedEdges.create(denoised.size(), denoised.type());
 
  
  Canny(denoised, detectedEdges, lowThresh, max_lowThresh,kernel);
  dst = Scalar::all(0);
  src.copyTo(dst, detectedEdges);  
  cv::namedWindow("Canny", CV_WINDOW_AUTOSIZE);
  cv::imshow("Canny", dst);

  vector< vector<cv::Point> > contours;
  vector<cv::Vec4i> hierarchy;
  cv::findContours(detectedEdges, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0,0));
  cv::Mat drawing = cv::Mat::zeros(detectedEdges.size(), CV_8UC3);
  for(int i = 0; i< contours.size(); i++)
  {
    cv::Scalar color = cv::Scalar( rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255) );
    drawContours( drawing, contours, i, color, 1, 8, vector<cv::Vec4i>(), 0, cv::Point() );
  }
  cv::namedWindow("contours", CV_WINDOW_AUTOSIZE);
  cv::imshow("contours", drawing);
*/
//  cv::namedWindow(srcName, CV_WINDOW_AUTOSIZE);
//  cv::imshow(srcName, src_gray);
//  cv::namedWindow(targetName, CV_WINDOW_AUTOSIZE);
//  cv::imshow(targetName, denoised);
  
  
//  cv::waitKey(0);
  return 0;

}
