if (${OpenCV_FOUND} AND ${Qt5Core_FOUND})
  add_executable(depthoffield
         ${GUI_SRCS}
         ${SUPERPIXELS_SRCS}
	 depthoffield.cc

  ) 
  target_link_libraries(depthoffield ${OpenCV_LIBS} ${GUI_LIBS})
endif()
