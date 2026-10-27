#pragma once
// Read-only adaptation of remollRunData v2 and remollTextFile v1. No remoll
// filesystem methods are needed. Fixed-width run counters avoid Windows long.
#include <TObject.h>
#include <TTimeStamp.h>
#include <cstring>
#include <string>
#include <vector>
struct filedata_t {
    char filename[255];
    char hashsum[255];
    TTimeStamp timestamp;
};
class remollTextFile:public TObject {
public:
    unsigned int fFilenameSize=0;
    char* fFilename=nullptr;
    unsigned long long int fBufferSize=0;
    char* fBuffer=nullptr;
    remollTextFile()=default;
    remollTextFile(const remollTextFile& other):TObject(other) {assign(other);}
    remollTextFile& operator=(const remollTextFile& other) {
        if(this!=&other){TObject::operator=(other);delete[] fFilename;delete[] fBuffer;assign(other);}return *this;
    }
    ~remollTextFile() override {delete[] fFilename;delete[] fBuffer;}
private:
    void assign(const remollTextFile& other) {
        fFilenameSize=other.fFilenameSize;fBufferSize=other.fBufferSize;
        fFilename=other.fFilename?new char[fFilenameSize]:nullptr;
        fBuffer=other.fBuffer?new char[fBufferSize]:nullptr;
        if(fFilename)std::memcpy(fFilename,other.fFilename,fFilenameSize);
        if(fBuffer)std::memcpy(fBuffer,other.fBuffer,fBufferSize);
    }
    ClassDefOverride(remollTextFile,1);
};
class remollRunData:public TObject {
public:
    std::vector<filedata_t> fMagData;
    remollTextFile fMacro;
    std::vector<remollTextFile> fGDMLFiles;
    Long64_t fNthrown=0;
    Long64_t fSeed=0;
    std::string fGitInfo;
    TTimeStamp fRunTime;
    std::string fRunPath;
    std::string fHostName;
    ClassDefOverride(remollRunData,3);
};
