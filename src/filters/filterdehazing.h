#pragma once

#include "filter.h"

class FilterDehazing : public Filter {
    static cv::Mat dehaze(const cv::Mat& image, const cv::Mat& transmittance, double min, double max)
    {
        cv::Mat colored_transmittance;
        cv::normalize(transmittance, colored_transmittance, min*255.0, max*255.0, cv::NORM_MINMAX, -1);
        cv::Mat1b atmosphere_mask  = (colored_transmittance <= 255.0*(min + 0.1*(max-min)));
        cv::Scalar atmosphere = cv::mean(image, atmosphere_mask);
        atmosphere/=255.0f;
        cv::cvtColor(colored_transmittance, colored_transmittance, CV_GRAY2RGB);
        cv::Mat num;
        cv::multiply(atmosphere, (colored_transmittance*(-1.0) + cv::Scalar(255.0,255.0,255.0)), num  , 1, CV_32FC3);
        cv::Mat imagef;
        image.convertTo(imagef, CV_32FC3);
        cv::subtract(imagef, num, num);
        cv::max(num, cv::Scalar(0.0,0.0,0.0), num);
        cv::Mat sol;
        cv::divide(num, colored_transmittance, sol, 255.0f, CV_8UC3);

        return sol;
    }

public:
    std::vector<PropagatedValue> propagatedValues() const override
    {
        return std::vector<PropagatedValue>{{
                PropagatedValue("Transmittance")
            }};

    }

    std::vector<FloatValue> floatValues() const override
    {
        return std::vector<FloatValue>{{
                FloatValue("Effect",0.0f,1.0f)
            }};
    }

    std::vector<Stroke> strokes() const override
    {
        return std::vector<Stroke>{{
                /*Stroke("Maximum fog", 0.25),
            Stroke("Medium fog",  0.50),
            Stroke("Minimum fog", 0.75),
            Stroke("No fog",      1.00),*/
            Stroke("maxmaximum fog",0.99),
            Stroke("maximum fog",0.75),
            Stroke("midmin fog",0.5),
            Stroke("minimum fog",0.25),
            Stroke("No fog",0.01),



            }};
    }



    cv::Mat apply(const cv::Mat& input_image,
                  const std::vector<std::shared_ptr<cv::Mat>>& propagated_values,
                  const std::vector<float>& float_values) const override
    {
        auto  transmittance = propagated_values[0];
        float intensity     = float_values[0];
        double min, max;
        cv::minMaxLoc(*transmittance, &min, &max);

        double min_t = (intensity*min) + (1.0 - intensity)*max;//0.05;
        double max_t = max;//0.95;
        
        return dehaze(input_image, *transmittance, min_t, max_t);
    }

};
