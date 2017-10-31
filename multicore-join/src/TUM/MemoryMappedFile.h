/*
 * MemoryMappedFile.h
 *
 *  Created on: Nov 8, 2012
 *      Author: Harald Lang
 */

#ifndef MEMORYMAPPEDFILE_H_
#define MEMORYMAPPEDFILE_H_

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

using namespace std;

class MemoryMappedFile {
private:
	void* map;
	char* filename;
	int fileDescriptor;
	uint64_t fileSize;

	MemoryMappedFile(const char* filename, int fileDescriptor, uint64_t fileSize)
				: fileDescriptor(fileDescriptor), fileSize(fileSize) {

		this->filename = new char[strlen(filename) + 1];
		strcpy(this->filename, filename);
		map = mmap(0, fileSize, PROT_READ | PROT_WRITE, MAP_SHARED, fileDescriptor, 0);
	};


public:
	~MemoryMappedFile() {
		closeFile();
		delete [] filename;
	}

	void* getMap() {
		return map;
	}

	void closeFile() {
		munmap(map, fileSize);
		close(fileDescriptor);
	}

	void deleteFile() {
		if (map != NULL) {
			closeFile();
		}
		remove(filename);
	}


	static MemoryMappedFile* createFile(string& filename, uint64_t fileSize) {
		// create a new file with the given size
		const char* fn = filename.c_str();
		int fd = open(fn, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
		lseek(fd, fileSize - 1, SEEK_SET);
		write(fd, "", 1);
		return new MemoryMappedFile(fn, fd, fileSize);
	}

	static MemoryMappedFile* openFile(string filename) {
		const char* fn = filename.c_str();
		int fd = open(fn, O_RDWR, (mode_t)0600);
		struct stat st;
		stat(fn, &st);
		uint64_t fileSize = st.st_size;
		return new MemoryMappedFile(fn, fd, fileSize);
	}

	template<typename T>
	static MemoryMappedFile* createFile(string& filename, uint64_t length)  {
		return createFile(filename, length * sizeof(T));
	}


};


#endif /* MEMORYMAPPEDFILE_H_ */
