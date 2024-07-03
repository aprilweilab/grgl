#ifndef GRGL_TESTING_UTILITIES_H
#define GRGL_TESTING_UTILITIES_H

#include <unistd.h>
#include <string>
#include <cstring>

static inline std::string writeTempFile(const std::string& contents, std::string ext = "") {
    constexpr size_t MAX_PATH = 256;
    char filename[MAX_PATH];
    strcpy(filename, "grgl_test_XXXXXX");
    int fd = mkstemp(filename);
	release_assert(fd > 0);
    release_assert(contents.size() == write(fd, contents.c_str(), contents.size()));
    fsync(fd);
    close(fd);

	if (!ext.empty()) {
		std::stringstream ss;
		ss << filename << ext;
		rename(filename, ss.str().c_str());
		return ss.str();
	}
    return std::string(filename);
}

static inline void remove_file(const std::string& filename) {
	unlink(filename.c_str());
}

#endif /* GRGL_TESTING_UTILITIES_H */