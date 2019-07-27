#pragma once

#include "filter.h"

class FilterDepthOfField : public Filter {

static std::vector<cv::Mat> create_blurred(const cv::Mat& input, double maxBlur, int nbins)
{
	std::vector<cv::Mat> blurred(nbins);
	blurred[0] = input;
	for (int i = 1; i<nbins;++i)
	{
		double blur = double(i)*maxBlur/double(nbins);
		int blur_size = int(blur);
		if (blur_size < 1) blur_size = 1;
		else if ((blur_size % 2) == 0) ++blur_size;
 		cv::GaussianBlur(input, blurred[i], cv::Size(blur_size,blur_size), 0);
	//	cv::imwrite("blurred_"+std::to_string(i)+".jpg",blurred[i]);
	}

	return blurred;
}


static cv::Mat interpolate_in_vector(const std::vector<cv::Mat>& source,
		const cv::Mat& interpolant)
{
    cv::Mat sol = cv::Mat::ones(source[0].rows, source[0].cols, source[0].type()); 
    cv::Mat_<cv::Vec3b> sol_  = sol;
    
    for( int i = 0; i < interpolant.rows; ++i) for( int j = 0; j < interpolant.cols; ++j )
    {
	float d = interpolant.at<float>(i,j);
	int   layer = trunc(d*float(source.size() - 1)); 
	if (layer >= (source.size()-1)) sol_(i,j) = source[source.size()-1].at<cv::Vec3b>(i,j);
	else {
		float t = std::min((d*float(source.size() - 1)) - float(layer), 1.0f);
		sol_(i,j) = source[layer  ].at<cv::Vec3b>(i,j)*(1.0f - t) +
		    	    source[layer+1].at<cv::Vec3b>(i,j)*t;
	}
    }

    sol = sol_;
    return sol;
}

static float blur_size_from_distance(float distance, float focal_distance, float focal_length, float aperture, bool linear)
{
	if (linear) return aperture*std::abs(distance - focal_distance)/focal_length;
	else return 
		aperture*(255.0 - focal_length)*std::abs(distance - focal_distance)/((255.0f - distance)*(focal_length - focal_distance));
}

static cv::Mat blur_image_depth(const cv::Mat& image, const cv::Mat& depth, 
		int nbins, float focal_distance, float focal_length, float aperture, bool linear) {
        
//	std::chrono::time_point<std::chrono::system_clock> start;
//        start = std::chrono::system_clock::now();

        cv::Mat sol;
        float ddepth = 1.0f/float(nbins);
	cv::Mat acc_mask;
	for (int i = 0; i<nbins;++i)
	{
		cv::Mat mask;
	        if (i == 0) mask = (depth >= -100.0f);
		else        mask = (depth <= (float(nbins - i)*ddepth));
		cv::cvtColor(mask, mask, CV_GRAY2RGB);
		cv::Mat masked;
		cv::bitwise_and(image, mask, masked);
		mask.convertTo(mask, CV_32F, 1.0f/255.0f);
		double blur = blur_size_from_distance((float(nbins - i - 1) + 0.5f)*ddepth,
				focal_distance, focal_length, aperture, linear);
		int blur_size = int(0.5*blur*float(image.cols));
		if (blur_size >= 1) {
			if ((blur_size % 2) == 0) ++blur_size;
 			cv::GaussianBlur(masked, masked, cv::Size(blur_size,blur_size), 0);
 			cv::GaussianBlur(mask,     mask, cv::Size(blur_size,blur_size), 0);
		}

		if (i == 0) {
			sol = masked.clone();
			acc_mask  = mask*(-1.0f) + cv::Scalar(1.0f,1.0f,1.0f); 
		} else {
			acc_mask  += mask*(-1.0f) + cv::Scalar(1.0f,1.0f,1.0f);
			cv::min(acc_mask, cv::Scalar(1.0f,1.0f,1.0f), acc_mask);
			cv::multiply(sol, acc_mask , sol, 1, CV_8U);
			cv::divide  (masked,      mask, masked        , 1, CV_8U);
			cv::multiply(masked, (acc_mask*(-1.0f) + cv::Scalar(1.0f,1.0f,1.0f)), masked   , 1, CV_8U);
		        sol += masked;	
		}
	}

   // std::cerr<<"Tiempo blur depth : "<<std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()<<std::endl;
    return sol;
}


//Blur image
static cv::Mat blur_image_focal_distance(const cv::Mat& image, const cv::Mat& depth, 
		int nbins, float focal_distance, float focal_length, float aperture, bool linear)
{
//    std::chrono::time_point<std::chrono::system_clock> start;

//    start = std::chrono::system_clock::now();
    auto blurred = create_blurred(image, aperture*255.0/focal_length, nbins);
//  imshow("blurred",blurred);
    cv::Mat distance_to_focus = cv::abs(depth - focal_distance);
    double min, max;
    cv::minMaxLoc(distance_to_focus, &min, &max);
    auto sol = interpolate_in_vector(blurred, distance_to_focus/max);  

    //std::cerr<<"Tiempo blur : "<<std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()<<std::endl;
    return sol;
    
}



public:
	std::vector<PropagatedValue> propagatedValues() const override
       	{    
		return std::vector<PropagatedValue>{{
			PropagatedValue("Depth")
		}};    
	
	}

	std::vector<FloatValue> floatValues() const override
       	{    
		return std::vector<FloatValue>{{
			FloatValue("Focal depth",0.0f,1.0f,0),
			FloatValue("Aperture",   0.0f,4.0f)
		}};    
	}


	std::vector<Stroke> strokes() const override
       	{    
		return std::vector<Stroke>{{
            Stroke("Furthest(1.0)", 1.0),
            Stroke("Far 2(0.67)",    0.67),
            Stroke("Far(0.34)",      0.34),
            Stroke("Close(0.15)",    0.15),
            Stroke("Near(0.01)",     0.01),
            //Stroke("user-defined", -1)
		}};
            /*return std::vector<Stroke>{{
                Stroke("Furthest", 0.01),
                Stroke("Far 2",    0.25),
                Stroke("Far",      0.44),
                Stroke("Close",    0.77),
                Stroke("Near",     1.00)
            }};*/
	}


	cv::Mat apply(const cv::Mat& input_image, 
			const std::vector<std::shared_ptr<cv::Mat>>& propagated_values,
			const std::vector<float>& float_values) const override
	{
		float focal_distance = float_values[0];
		float aperture       = float_values[1];
        	auto  depth          = (1-*propagated_values[0]);

		float focal_length = focal_distance + 20;
        
	        return blur_image_depth(input_image, depth, 6, focal_distance, focal_length, aperture, true);
	}

};
