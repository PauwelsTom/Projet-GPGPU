all:
	cmake -S . -B build --preset release -D USE_CUDA=ON
	@echo -e "\n\n\n"
	cmake --build build
	@echo -e "\n\n\n"
	wget https://gstreamer.freedesktop.org/media/sintel_trailer-480p.webm
	@echo -e "\n\n\n"
	export GST_PLUGIN_PATH=$(shell pwd)
	@echo -e "\n\n\n"
	touch libgstcudafilter.so 
	rm libgstcudafilter.so
	@echo -e "\n\n\n"
	ln -s ./build/libgstcudafilter-cpp.so libgstcudafilter.so
	@echo -e "\n\n\n"
	gst-launch-1.0 uridecodebin uri=file://$(shell pwd)/sintel_trailer-480p.webm ! videoconvert ! "video/x-raw, format=(string)RGB" ! cudafilter ! videoconvert ! video/x-raw, format=I420 ! x264enc ! mp4mux ! filesink location=video.mp4


run:
	cmake --build build
	@echo -e "\n\n\n"
	gst-launch-1.0 uridecodebin uri=file://$(shell pwd)/sintel_trailer-480p.webm ! videoconvert ! "video/x-raw, format=(string)RGB" ! cudafilter ! videoconvert ! video/x-raw, format=I420 ! x264enc ! mp4mux ! filesink location=video.mp4
