#pragma once

#include <vector>
#include <list>
#include <string>
#include <tuple>
#include <opencv2/core/core.hpp>


class Filter {
protected:
    static std::string type2str(int type)
    {
		std::string r;

		  unsigned char depth = type & CV_MAT_DEPTH_MASK;
		  unsigned char chans = 1 + (type >> CV_CN_SHIFT);

		  switch ( depth ) {
		    case CV_8U:  r = "8U"; break;
		    case CV_8S:  r = "8S"; break;
		    case CV_16U: r = "16U"; break;
		    case CV_16S: r = "16S"; break;
		    case CV_32S: r = "32S"; break;
		    case CV_32F: r = "32F"; break;
		    case CV_64F: r = "64F"; break;
		    default:     r = "User"; break;
		  }

		  r += "C";
		  r += (chans+'0');

		  return r;
	}

	static void print_info(const char* name, const cv::Mat& image) 
	{
		double min, max;
		cv::minMaxLoc(image, &min, &max);

		std::cerr<<std::setw(20)<<name<<" ("<<type2str(image.type())<<")\t- "<<image.cols<<"x"<<image.rows<<"x"<<image.channels()<<" - ["<<min<<","<<max<<"]"<<std::endl;
	}

public:
    class PropagatedValue
    {
		std::string _name;
		float _default_value;
	public:	
		PropagatedValue(const std::string& name = std::string(), float default_value = 0.5f) :
			_name(name), _default_value(default_value) { }

		PropagatedValue(const char* name, float default_value = 0.5f) :
			PropagatedValue(std::string(name), default_value) { }

		const std::string& name() const { return _name; }
		float defaultValue() const { return _default_value; }
	};

    class FloatValue
    {
		std::string _name;
		float _min, _max;
		int _pickable_from_channel;
	public:
		FloatValue(const std::string& name = std::string(), float min = 0.0f, float max = 1.0f, int pickable_from_channel = -1) :
			_name(name), _min(min), _max(max), _pickable_from_channel(pickable_from_channel) { }

		FloatValue(const char* name, float min, float max, int pickable_from_channel = -1) : 
			FloatValue(std::string(name), min, max, pickable_from_channel) { }
		const std::string& name() const { return _name; }
		float min() const { return _min; }
		float max() const { return _max; }
		bool isPickable() const { return _pickable_from_channel >= 0; }
		int  channelToPickFrom() const { return _pickable_from_channel; }
	};

    class Stroke
    {
		std::string _name;
		std::list<std::tuple<float, int>> _values;

	public:
		Stroke(const std::string& name, const std::list<std::tuple<float, int>>& values) :
			_name(name), _values(values) { }

		Stroke(const std::string& name = std::string(), float value=0.5, int channel=0) :
			Stroke(name, std::list<std::tuple<float, int>>{{std::make_tuple(value, channel)}}) { }

        Stroke(const std::string& name, float value1, int channel1, float value2, int channel2)
            :
            Stroke(name, std::list<std::tuple<float, int>>
            {
                {
                    std::make_tuple(value1, channel1), std::make_tuple(value2, channel2)
                }
            }
            ) { }
		Stroke(const char* name, const std::list<std::tuple<float, int>>& values) :
			Stroke(std::string(name), values) { }
		Stroke(const char* name, float value, int channel=0) :
			Stroke(std::string(name), value, channel) { }
		Stroke(const char* name, float value1, int channel1, float value2, int channel2) :
			Stroke(std::string(name), value1, channel1, value2, channel2) { }

		const std::string& name() const { return _name; }
		const std::list<std::tuple<float, int>>& values() const { return _values; }
	};

	virtual std::vector<PropagatedValue> propagatedValues() const
       	{    return std::vector<PropagatedValue>(0);    }
	virtual std::vector<FloatValue> floatValues() const
       	{    return std::vector<FloatValue>(0);    }
	virtual std::vector<Stroke> strokes() const
       	{    return std::vector<Stroke>(0);    }

	virtual cv::Mat apply(const cv::Mat& input_image, 
			const std::vector<std::shared_ptr<cv::Mat>>& propagated_values,
			const std::vector<float>& float_values) const = 0;

};
