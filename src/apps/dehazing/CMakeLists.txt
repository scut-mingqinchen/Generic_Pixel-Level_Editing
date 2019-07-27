if (${OpenCV_FOUND} AND ${Qt5Core_FOUND})
  add_executable(dehazing
         ${GUI_SRCS}
         ${SUPERPIXELS_SRCS}
	 dehazing.cc
  ) 
  target_link_libraries(dehazing ${OpenCV_LIBS} ${GUI_LIBS})
endif()
