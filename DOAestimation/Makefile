CXX = g++
CXXFLAGS = -Wall -std=c++11 -g -O0
LDFLAGS = -lmatio -lfftw3

# Common source files used by all executables
COMMON_SRCS = DoasCurrentFrame.cpp \
              locate_zones.cpp

# Block v4 specific files
BLOCK_V4_SPECIFIC = DOA_Histogram.cpp \
                    processHist.cpp

# Block v3 specific files
BLOCK_V3_SPECIFIC = DOA_Histogram.cpp

# Block v4 executable
BLOCK_V4_SRCS = run_DRACOSS_1D_Blockbased_v4.cpp $(COMMON_SRCS) $(BLOCK_V4_SPECIFIC)
BLOCK_V4_OBJS = $(BLOCK_V4_SRCS:.cpp=.o)
BLOCK_V4_TARGET = doa_block_v4

# Block v3 executable
BLOCK_V3_SRCS = run_DRACOSS_1D_Allframes_Blockbased_script_v03.cpp $(COMMON_SRCS) $(BLOCK_V3_SPECIFIC)
BLOCK_V3_OBJS = $(BLOCK_V3_SRCS:.cpp=.o)
BLOCK_V3_TARGET = doa_block_v3

# Single frame executable
SINGLE_FRAME_SRCS = run_DRACOSS_1D_singleframe_script.cpp $(COMMON_SRCS)
SINGLE_FRAME_OBJS = $(SINGLE_FRAME_SRCS:.cpp=.o)
SINGLE_FRAME_TARGET = doa_single_frame

# All frames executable
ALL_FRAMES_SRCS = run_DRACOSS_1D_Allframes.cpp $(COMMON_SRCS)
ALL_FRAMES_OBJS = $(ALL_FRAMES_SRCS:.cpp=.o)
ALL_FRAMES_TARGET = doa_all_frames

# Default target
all: $(BLOCK_V4_TARGET) $(BLOCK_V3_TARGET) $(SINGLE_FRAME_TARGET) $(ALL_FRAMES_TARGET)

# Block v4
$(BLOCK_V4_TARGET): $(BLOCK_V4_OBJS)
	$(CXX) $(BLOCK_V4_OBJS) -o $(BLOCK_V4_TARGET) $(LDFLAGS)

# Block v3
$(BLOCK_V3_TARGET): $(BLOCK_V3_OBJS)
	$(CXX) $(BLOCK_V3_OBJS) -o $(BLOCK_V3_TARGET) $(LDFLAGS)

# Single frame
$(SINGLE_FRAME_TARGET): $(SINGLE_FRAME_OBJS)
	$(CXX) $(SINGLE_FRAME_OBJS) -o $(SINGLE_FRAME_TARGET) $(LDFLAGS)

# All frames
$(ALL_FRAMES_TARGET): $(ALL_FRAMES_OBJS)
	$(CXX) $(ALL_FRAMES_OBJS) -o $(ALL_FRAMES_TARGET) $(LDFLAGS)

# Pattern rule for object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean all
clean:
	rm -f *.o $(BLOCK_V4_TARGET) $(BLOCK_V3_TARGET) $(SINGLE_FRAME_TARGET) $(ALL_FRAMES_TARGET)

# Run targets
run_block_v4: $(BLOCK_V4_TARGET)
	./$(BLOCK_V4_TARGET)

run_block_v3: $(BLOCK_V3_TARGET)
	./$(BLOCK_V3_TARGET)

run_single_frame: $(SINGLE_FRAME_TARGET)
	./$(SINGLE_FRAME_TARGET)

run_all_frames: $(ALL_FRAMES_TARGET)
	./$(ALL_FRAMES_TARGET)

.PHONY: all clean run_block_v4 run_block_v3 run_single_frame run_all_frames 