#pragma once
#define CHUNK_HEIGH 120
#include <string>
#include <vector>
#include <iostream>
#include "gdal_priv.h"
using namespace std;

class VI
{
public:
	VI();
	~VI();
	void Compute_NDVI(vector<string> fileName, char * OutFileName, char* MCD12Q1);
	bool getTimeChar(char timeChar[]);
	

private:
	void readModisPara(char* modis_path, short* modis_Para, int bandName);

	bool getGeo(char* modis_path, string & wkt, string & adfGeoTransform);

	void VI_Composite_1000m_NDVI(float* solarZenithAngle, float* viewZenithAngle, float* RAAngle, float* RedData,
		float* NirData, float* NDVIData, unsigned char* cloudMask,
		unsigned char landCover, short &new_ndvi, short &new_ndviqc, int size);

	vector<string> getNumFile(vector<string> fileName, int *senorID);

	bool fileDetection(char* filePath);
	/**
	*@brief 读取数据
	*/
	bool readFileData(vector<string> fileName, int *senorID, int offx, int offy, int pafsizeX, int pafsizeY);
	/**
	*@brief 读取MODIS数据
	*/
	bool readMODISData(char* fileName, int i, int offx, int offy, int pafsizeX, int pafsizeY);
	/**
	*@brief 读取FY_VIRR数据
	*/
	bool readFYVIRRData(char* fileName, int i, int offx, int offy, int pafsizeX, int pafsizeY);
	/**
	*@brief 读取FY_MERSI数据
	*/
	bool readFYMERSIData(char* fileName, int i, int offx, int offy, int pafsizeX, int pafsizeY);

	bool getMemory(int fileNum, int pafsizeX, int pafsizeY);

	bool freeMemory(int fileNum);
	/**
	*@brief 获得网格编号
	*/
	string getGridNum(char * OutFileName);
	/**
	*@brief 获得分类数据
	*/
	string getLandCoverData(string GridNumStr, string classPath);
	/**
	*@brief 获得分类数据的路径
	*/
//	bool getTmpPath(string &classPath);
	/**
	*@brief 获得所有数据的路径
	*/
	vector<string> getMCDPath(const char* ss);
	/**
	*@brief 计算均值
	*/
	float mean(vector<float> data);
	/**
	*@brief 计算方差
	*/
	float std(vector<float> data);
	int getMax(vector<float> data);

	void readLandCover(char* LandCoverDataPath, unsigned char* LandCoverData);
	/**
	*@brief 从完整路径中得到文件名，不带扩展名
	*@param  strFullPath    	[in] 完整路径
	*@param  result  	    [in] 文件名
	*@return
	*/
	const char* getFilaNameFromPath(const char* strFullPath);
	vector<string> getMCDPath(const char* ss, const char* splitC);

private:
	int nXSize;
	int nYSize;
	unsigned short** BlueData;
	unsigned short** RedData;
	unsigned short** NirData;
	short** solarZenithData;
	short** viewZenithData;
	short** solarAzimuthData;
	short** viewAzimuthData;
	unsigned char** cloudData;
	
};

