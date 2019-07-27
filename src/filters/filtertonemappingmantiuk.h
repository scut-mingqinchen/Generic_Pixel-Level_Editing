#pragma once

#include "filter.h"

class FilterTonemappingMantiuk : public Filter {

static float luminance_curve(float luminance, float b, float c, float dl, float dh)
{
	if (luminance<=0.0f) return 0.0f;
	float ll = std::log(luminance);
	float al = (c*dl - 1.0)/dl;
	float ah = (c*dh - 1.0)/dh;
	if (ll<(b-dl))      return 0.0;
	else if (ll<b)      return 0.5*c*(ll - b)/(1 - al*(ll-b)) + 0.5;
	else if (ll<(b+dh)) return 0.5*c*(ll - b)/(1 + ah*(ll-b)) + 0.5;
	else                return 1.0;
}

public:
//Trying this: https://www.cl.cam.ac.uk/~rkm38/pdfs/mantiuk08mgtmo.pdf
static cv::Mat tonemap(const cv::Mat& image, const cv::Mat& brightness, const cv::Mat& contrast, float dl, float dh)
{
    cv::Mat sol;

    std::vector<cv::Mat> channels;
    cv::split(image, channels);
    cv::Mat& B = channels[0];
    cv::Mat& G = channels[1];
    cv::Mat& R = channels[2];
    cv::Mat luminance = 0.27*R + 0.67*G + 0.06*B;

    //These below are for the color correction.
    cv::Mat r =  R/luminance;
    cv::Mat g =  G/luminance;
    cv::Mat b =  B/luminance;

    cv::Mat corrected_luminance(luminance.rows, luminance.cols, CV_32F);
    for(int i=0; i<luminance.rows; i++) 
    	for(int j=0; j<luminance.cols; j++) 
		corrected_luminance.at<float>(i,j) = luminance_curve(luminance.at<float>(i,j), brightness.at<float>(i,j), contrast.at<float>(i,j), dl, dh);
    
    cv::multiply(r, corrected_luminance, R, 255.0, CV_32F); 
    cv::multiply(g, corrected_luminance, G, 255.0, CV_32F); 
    cv::multiply(b, corrected_luminance, B, 255.0, CV_32F); 
    cv::merge(channels, sol);
    
    cv::max(sol, cv::Scalar(0.0,0.0,0.0), sol);
    cv::min(sol, cv::Scalar(255.0,255.0,255.0), sol);
    sol.convertTo(sol, CV_8UC3);
    return sol;
}


public:
	std::vector<PropagatedValue> propagatedValues() const override
       	{    
		return std::vector<PropagatedValue>{{
			PropagatedValue("Local brightness"),
			PropagatedValue("Local contrast")
		}};    
	
	}

	std::vector<FloatValue> floatValues() const override
       	{    
		return std::vector<FloatValue>{{
			FloatValue("Brightness",  -2.0f,4.0f),
			FloatValue("Contrast",     0.0f,2.0f),
			FloatValue("Local brightness", 0.0f,4.0f),
			FloatValue("Local contrast",   0.0f,4.0f)
		}};    
	}

	std::vector<Stroke> strokes() const override
       	{    
		return std::vector<Stroke>{{
			Stroke("Preserve",       0.50,0, 0.50,1),
			Stroke("Darken",         0.10,0),
			Stroke("Light up",       1.00,0),
			Stroke("Reduce contrast",   0.05,1),
			Stroke("Increase contrast", 1.00,1)
		}};

	}



	cv::Mat apply(const cv::Mat& input_image, 
			const std::vector<std::shared_ptr<cv::Mat>>& propagated_values,
			const std::vector<float>& float_values) const override
	{	
		float brightness              = float_values[0];
		float contrast                = float_values[1];
		float effect_local_brightness = float_values[2];
		float effect_local_contrast   = float_values[3];
		auto  local_brightness        = propagated_values[0];
		auto  local_contrast          = propagated_values[1];

		cv::Mat b;
            cv::exp(effect_local_brightness*(cv::Scalar(0.5)-*local_brightness), b);
        b*=-brightness;
        
		cv::Mat c;
            cv::exp(effect_local_contrast*(*local_contrast - cv::Scalar(0.5)), c);
        c*=contrast;

		double min, max;
		cv::minMaxLoc(input_image, &min, &max);
		float logmin = std::log(std::max(min, 1.e-2));
		float logmax = std::log(max);
        
        return tonemap(input_image, b, c, 2, 1);
	}

};
