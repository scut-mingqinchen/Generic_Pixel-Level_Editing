if (${OpenCV_FOUND} AND ${Qt5Core_FOUND})
  add_executable(tonemapping
         ${GUI_SRCS}
         ${SUPERPIXELS_SRCS}
	 tonemapping.cc
  ) 
  target_link_libraries(tonemapping ${OpenCV_LIBS} ${GUI_LIBS})
endif()
