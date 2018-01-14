#include "VI.h"
#include "maparams.h"
#include <fstream>
#include "QualityCal.h"

VI::VI()
{
}


VI::~VI()
{
}

void VI::Compute_NDVI(vector <string>  fileName, char * OutFileName, char * MCD12Q1)
{
	GDALAllRegister();
	string GridNumStr = getGridNum(OutFileName);   //获取条带号：H??V??
//	vector<string> fileName = getMCDPath(strFileName);
	int* senorID = new int[3];//各个传感器有的数据数；
	vector<string> filenameUse = getNumFile(fileName, senorID);//获得各传感器有效的文件数据，并且排序
	int fileNum = filenameUse.size();
	string projection;
	string aodGeoTransform;
	getGeo(const_cast<char*>(filenameUse[0].c_str()), projection, aodGeoTransform);
	/////create output file??
	string strFileNameOnly;   //用作写入元数据
	for (int i = 0; i < fileNum; i++)
	{
		string NameOnly(getFilaNameFromPath(filenameUse[i].c_str()));
		cout << filenameUse[i] << endl;
		strFileNameOnly += NameOnly;
		strFileNameOnly += ";";
	}
	char time[20];
	getTimeChar(time);

	nXSize = 1200;
	nYSize = 1200;
	maparams ss;
	ss.Red = new float[fileNum];   //composite numbers
	ss.NIR = new float[fileNum];
	ss.Blue = new float[fileNum];
	ss.NDVI = new float[fileNum];
	ss.EVI = new float[fileNum];
	ss.viewZenith = new float[fileNum];
	ss.solarZenith = new float[fileNum];
	ss.RaAngle = new float[fileNum];
	ss.cloudMask = new unsigned char[fileNum];

//	int* qcref_new = new int[fileNum];  

	short* m_NewNdvi = new short[nXSize*CHUNK_HEIGH];
	short* m_NewEvi = new short[nXSize*CHUNK_HEIGH];
	short* m_QANdvi = new short[nXSize*CHUNK_HEIGH];
	short* m_QAEvi = new short[nXSize*CHUNK_HEIGH];
	//short* modisPara1 = new short[1200 * 1200 * 3];
	//short* modisPara2 = new short[1200 * 1200 * 3];
	unsigned char* LandCoverData = new unsigned char[nXSize*nYSize];
	readLandCover(MCD12Q1, LandCoverData);
	int offx = 0, offy = 0;
	int pafsizeX = nXSize, pafsizeY = CHUNK_HEIGH;
	int num = (nYSize - 1) / CHUNK_HEIGH + 1;     //subset num
	getMemory(fileNum, pafsizeX, pafsizeY);

	for (int kk = 0; kk < num; kk++)
	{
		cout << kk << endl;
		offy = kk*CHUNK_HEIGH;
		if (kk == (num - 1))
		{
			pafsizeY = (nYSize - 1) % CHUNK_HEIGH + 1;
		}
		readFileData(filenameUse, senorID, offx, offy, pafsizeX, pafsizeY);
		for (int i = 0; i < pafsizeY; i++)
		{
			for (int j = 0; j < pafsizeX; j++)
			{
				int day = fileNum;
				int countEffective = 0;
				for (int dd = 0; dd < day; dd++)
				{
					if (BlueData[dd][i*pafsizeX + j] == 65535 || RedData[dd][i*pafsizeX + j] == 65535 ||
						NirData[dd][i*pafsizeX + j] == 65535 || solarZenithData[dd][i*pafsizeX + j] == -32767 ||
						viewZenithData[dd][i*pafsizeX + j] == -32767 || viewAzimuthData[dd][i*pafsizeX + j] == -32767 ||
						solarAzimuthData[dd][i*pafsizeX + j] == -32767 || cloudData[dd][i*pafsizeX + j] == 255)
					{
						ss.Blue[dd] = 0;
						ss.Red[dd] = 0;
						ss.NIR[dd] = 0;
						ss.solarZenith[dd] = 0;
						ss.viewZenith[dd] = 0;
						ss.RaAngle[dd] = 0;
						ss.cloudMask[dd] = 0;
					}
					else
					{
						ss.Blue[dd] = (BlueData[dd][i*pafsizeX + j]) / 10000.0;
						ss.Red[dd] = (RedData[dd][i*pafsizeX + j]) / 10000.0;
						ss.NIR[dd] = (NirData[dd][i*pafsizeX + j]) / 10000.0;
						ss.solarZenith[dd] = (solarZenithData[dd][i*pafsizeX + j]) / 100.0;
						ss.viewZenith[dd] = (viewZenithData[dd][i*pafsizeX + j]) / 100.0;
						ss.RaAngle[dd] = (viewAzimuthData[dd][i*pafsizeX + j] - solarAzimuthData[dd][i*pafsizeX + j]) / 100.0;
						ss.cloudMask[dd] = cloudData[dd][i*pafsizeX + j];

					}
					if (ss.Red[dd] <= 0 || ss.Red[dd] >= 1 || ss.NIR[dd] <= 0 || ss.NIR[dd] >= 1)
					{
						ss.NDVI[dd] = 0;
						ss.EVI[dd] = 0;
					}
					else
					{
						countEffective++;
						float ndvitmp = (ss.NIR[dd] - ss.Red[dd]) / (ss.NIR[dd] + ss.Red[dd]);
						float evitmp = G1*(ss.NIR[dd] - ss.Red[dd]) / (ss.NIR[dd] + C1*ss.Red[dd] - C2*ss.Blue[dd] + L1);
						if (ndvitmp<0 || ndvitmp>1)
						{
							ss.NDVI[dd] = 0;
						}
						else
						{
							ss.NDVI[dd] = ndvitmp;
						}
						if (evitmp<0 || evitmp>1)
						{
							ss.EVI[dd] = 0;
						}
						else
						{
							ss.EVI[dd] = evitmp;
						}
					}//else
				}//file iter
				if (countEffective != 0)
				{
			//		memset(qcref_new, 5, sizeof(int)*day);
			//		qc.getQcref(ss.Red, ss.NIR, ss.solarZenith, ss.viewZenith, ss.RaAngle, MCD43Data, LandCoverData[(i + offy)*nXSize + j + offx],
			//			qcref_new, senorID[0], day, ss.cloudMask);

					short new_ndvi = 0, new_ndviqc = 0;
					short new_evi = 0, new_eviqc = 0;
					VI_Composite_1000m_NDVI(ss.solarZenith, ss.viewZenith, ss.RaAngle, ss.Red, ss.NIR, ss.NDVI, ss.cloudMask,
						 LandCoverData[(i + offy)*nXSize + j + offx], new_ndvi, new_ndviqc, day);
			//		VI_Composite_1000m_EVI(ss.solarZenith, ss.viewZenith, ss.RaAngle, ss.Blue, ss.Red, ss.NIR, ss.EVI, ss.cloudMask,
			//			qcref_new, LandCoverData[(i + offy)*nXSize + j + offx], new_evi, new_eviqc, day);
					m_NewNdvi[i*pafsizeX + j] = new_ndvi;
					m_QANdvi[i*pafsizeX + j] = new_ndviqc;
					m_NewEvi[i*pafsizeX + j] = new_evi;
					m_QAEvi[i*pafsizeX + j] = new_eviqc;
				}
				else
				{
					m_NewNdvi[i*pafsizeX + j] = -32767;
					m_QANdvi[i*pafsizeX + j] = -32767;
					m_NewEvi[i*pafsizeX + j] = -32767;
					m_QAEvi[i*pafsizeX + j] = -32767;
				}
			}
		}

		/*pDataSet_NDVI->WriteRaster(offx, offy, pafsizeX, pafsizeY, m_NewNdvi, DFAL_DT_Int16, 0, NULL);
		pDataSet_NDVIqc->WriteRaster(offx, offy, pafsizeX, pafsizeY, m_QANdvi, DFAL_DT_Int16, 0, NULL);
		pDataSet_EVI->WriteRaster(offx, offy, pafsizeX, pafsizeY, m_NewEvi, DFAL_DT_Int16, 0, NULL);
		pDataSet_EVIqc->WriteRaster(offx, offy, pafsizeX, pafsizeY, m_QAEvi, DFAL_DT_Int16, 0, NULL);*/

	}

}

bool VI::getTimeChar(char timeChar[])
{
	time_t now; //实例化time_t结构
	struct tm *timenow; //实例化tm结构指针
	time(&now);
	timenow = localtime(&now);
	char year[5], day[4], hour[3], min[3], sec[3];
	sprintf(year, "%d", timenow->tm_year + 1900);
	strcpy(timeChar, year);

	//判断天数*******************************
	if (timenow->tm_yday<10)
	{
		sprintf(day, "%d", timenow->tm_yday);
		strcat(timeChar, "00");
		strcat(timeChar, day);
	}
	else if (timenow->tm_yday<100)
	{
		sprintf(day, "%d", timenow->tm_yday);
		strcat(timeChar, "0");
		strcat(timeChar, day);
	}
	else
	{
		sprintf(day, "%d", timenow->tm_yday);
		strcat(timeChar, day);
	}
	//判断小时*******************************
	if (timenow->tm_hour<10)
	{
		sprintf(hour, "%d", timenow->tm_hour);
		strcat(timeChar, "0");
		strcat(timeChar, hour);
	}
	else
	{
		sprintf(hour, "%d", timenow->tm_hour);
		strcat(timeChar, hour);
	}
	//判断分钟*******************************
	if (timenow->tm_min<10)
	{
		sprintf(min, "%d", timenow->tm_min);
		strcat(timeChar, "0");
		strcat(timeChar, min);
	}
	else
	{
		sprintf(min, "%d", timenow->tm_min);
		strcat(timeChar, min);
	}
	//判断秒********************************
	if (timenow->tm_sec<10)
	{
		sprintf(sec, "%d", timenow->tm_sec);
		strcat(timeChar, "0");
		strcat(timeChar, sec);
	}
	else
	{
		sprintf(sec, "%d", timenow->tm_sec);
		strcat(timeChar, sec);
	}

	return true;
	
}




void VI::readModisPara(char* modis_path, short* modis_Para, int bandName)
{

}


bool VI::getGeo(char * modis_path, string & wkt, string & adfGeoTransform)
{
	GDALDataset* podataset = (GDALDataset*)GDALOpen(modis_path, GA_ReadOnly);
	const char* adf = GDALGetMetadataItem((GDALDatasetH)podataset, "ProjectionPara", NULL);
	const char * proj1= GDALGetMetadataItem((GDALDatasetH)podataset, "ProjectionStr", NULL);
	wkt = (string)proj1;
	adfGeoTransform = (string)adf;

//	const char* test = GDALGetMetadataItem((GDALDatasetH)podataset, "SUBDATASET_1_NAME" ,"Subdatasets");
	return true;
}

void VI::VI_Composite_1000m_NDVI(float * solarZenithAngle, float * viewZenithAngle, float * RAAngle, float * RedData, 
	float * NirData, float * NDVIData, unsigned char * cloudMask, unsigned char landCover, short & new_ndvi, short & new_ndviqc, int size)
{
	if (landCover == 0)
	{
		new_ndvi = 0;
		new_ndviqc = 128;
	}
	else
	{
		float tn_ndvi, tn_ndviqc;
		QualityCal qcndvi(RedData,NirData,solarZenithAngle, viewZenithAngle,RAAngle,landCover,0,size,cloudMask);
		qcndvi.Process(tn_ndvi,tn_ndviqc);
	}
}

vector<string> VI::getNumFile(vector<string> fileName, int * senorID)
{
	vector<string> fileTname;
	for (int i = 0; i < 3; i++)
	{
		senorID[i] = 0;
	}
	for (int i = 0; i < fileName.size(); i++)
	{
		string namePath = (string)fileName[i];
		string name = namePath.substr(namePath.find_last_of("/") + 1);
		if ((name.find("MODIS") != string::npos || name.find("MOD021KM") != string::npos || name.find("MYD021KM") != string::npos) && (fileDetection(const_cast<char*>(fileName[i].c_str())) == true))
		{
			fileTname.push_back(fileName[i]);
			senorID[0]++;
		}

	}
	for (int i = 0; i < fileName.size(); i++)
	{
		string namePath = (string)fileName[i];
		string name = namePath.substr(namePath.find_last_of("/") + 1);
		if (name.find("MERSI") != string::npos && (fileDetection(const_cast<char*>(fileName[i].c_str())) == true))
		{
			fileTname.push_back(fileName[i]);
			senorID[1]++;
		}
	}
	for (int i = 0; i < fileName.size(); i++)
	{
		string namePath = (string)fileName[i];
		string name = namePath.substr(namePath.find_last_of("/") + 1);
		if (name.find("VIRR") != string::npos && (fileDetection(const_cast<char*>(fileName[i].c_str())) == true))
		{
			fileTname.push_back(fileName[i]);
			senorID[2]++;
		}
	}
	return fileTname;
}

bool VI::fileDetection(char * filePath)
{
	
	GDALDataset* podataset = (GDALDataset*)GDALOpen(filePath, GA_ReadOnly);
	if (podataset == NULL)
	{
		return false;
	}
	char** sublist = GDALGetMetadata((GDALDatasetH)podataset, "SUBDATASETS");
	string subDataset1Name = sublist[0];
	subDataset1Name = subDataset1Name.substr(subDataset1Name.find_first_of("=") + 1);
	GDALDataset* subDataset = (GDALDataset*)GDALOpen(subDataset1Name.c_str(), GA_ReadOnly);//打开该数据  
	if (subDataset == NULL)
	{
		return false;
	}
	GDALClose((GDALDatasetH)subDataset);
	GDALClose((GDALDatasetH)podataset);
	return true;
}

bool VI::readFileData(vector<string> fileName, int * senorID, int offx, int offy, int pafsizeX, int pafsizeY)
{
	for (int i = 0; i < fileName.size(); i++)
	{
		if (i < senorID[0])//MOD数据,MYD数据
		{
			readMODISData(const_cast<char*>(fileName[i].c_str()), i, offx, offy, pafsizeX, pafsizeY);
		}
		else if (i < (senorID[0] + senorID[1])) //FY-MERSI
		{
			readFYMERSIData(const_cast<char*>(fileName[i].c_str()), i, offx, offy, pafsizeX, pafsizeY);
		}
		else if (i < (senorID[0] + senorID[1] + senorID[2])) //FY-VIRR
		{
			readFYVIRRData(const_cast<char*>(fileName[i].c_str()), i, offx, offy, pafsizeX, pafsizeY);
		}
	}
	return true;
}

bool VI::readMODISData(char * fileName, int i, int offx, int offy, int pafsizeX, int pafsizeY)
{
	GDALDataset* pdataset = (GDALDataset*)GDALOpen(fileName, GA_ReadOnly);
	if (pdataset == NULL)
	{
		cout << "文件打开失败:" << fileName << endl;
		return false;//文件打开失败
	}

	GDALDataset* pDataSetSA = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_1_NAME", "Subdatasets"),GA_ReadOnly);
	GDALDataset* pDataSetSZ = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_2_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetVA = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_3_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetVZ = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_4_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetBlue = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_7_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetRed = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_5_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetNIR = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_6_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetCM = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_12_NAME", "Subdatasets"), GA_ReadOnly);
	if (pDataSetSA == NULL || pDataSetSZ == NULL || pDataSetVA == NULL || pDataSetVZ == NULL ||
		pDataSetRed == NULL || pDataSetNIR == NULL || pDataSetCM == NULL)
	{
		cout << "文件数据集打开失败:" << fileName << endl;
		return false;//文件数据集打开失败
	}
	pDataSetSA->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, solarAzimuthData[i],
		pafsizeX, pafsizeY, pDataSetSA->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetSZ->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, solarZenithData[i],
		pafsizeX, pafsizeY, pDataSetSZ->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetVA->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, viewAzimuthData[i],
		pafsizeX, pafsizeY, pDataSetVA->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetVZ->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, viewZenithData[i],
		pafsizeX, pafsizeY, pDataSetVZ->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetBlue->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, BlueData[i],
		pafsizeX, pafsizeY, pDataSetBlue->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetRed->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, RedData[i],
		pafsizeX, pafsizeY, pDataSetRed->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetNIR->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, NirData[i],
		pafsizeX, pafsizeY, pDataSetNIR->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetCM->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, cloudData[i],
		pafsizeX, pafsizeY, pDataSetCM->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);

	GDALClose((GDALDatasetH)pDataSetSA);
	GDALClose((GDALDatasetH)pDataSetSZ);
	GDALClose((GDALDatasetH)pDataSetVA);
	GDALClose((GDALDatasetH)pDataSetVZ);
	GDALClose((GDALDatasetH)pDataSetBlue);
	GDALClose((GDALDatasetH)pDataSetRed);
	GDALClose((GDALDatasetH)pDataSetNIR);
	GDALClose((GDALDatasetH)pdataset);
	return true;
}

bool VI::readFYVIRRData(char * fileName, int i, int offx, int offy, int pafsizeX, int pafsizeY)
{
	return false;
}

bool VI::readFYMERSIData(char * fileName, int i, int offx, int offy, int pafsizeX, int pafsizeY)
{
	GDALDataset* pdataset = (GDALDataset*)GDALOpen(fileName, GA_ReadOnly);
	if (pdataset == NULL)
	{
		cout << "文件打开失败:" << fileName << endl;
		return false;//文件打开失败
	}

	GDALDataset* pDataSetSA = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_1_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetSZ = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_2_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetVA = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_3_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetVZ = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_4_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetBlue = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_5_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetRed = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_7_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetNIR = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_8_NAME", "Subdatasets"), GA_ReadOnly);
	GDALDataset* pDataSetCM = (GDALDataset*)GDALOpen(pdataset->GetMetadataItem("SUBDATASET_10_NAME", "Subdatasets"), GA_ReadOnly);
	if (pDataSetSA == NULL || pDataSetSZ == NULL || pDataSetVA == NULL || pDataSetVZ == NULL ||
		pDataSetRed == NULL || pDataSetNIR == NULL || pDataSetCM == NULL)
	{
		cout << "文件数据集打开失败:" << fileName << endl;
		return false;//文件数据集打开失败
	}
	pDataSetSA->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, solarAzimuthData[i],
		pafsizeX, pafsizeY, pDataSetSA->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetSZ->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, solarZenithData[i],
		pafsizeX, pafsizeY, pDataSetSZ->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetVA->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, viewAzimuthData[i],
		pafsizeX, pafsizeY, pDataSetVA->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetVZ->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, viewZenithData[i],
		pafsizeX, pafsizeY, pDataSetVZ->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetBlue->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, BlueData[i],
		pafsizeX, pafsizeY, pDataSetBlue->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetRed->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, RedData[i],
		pafsizeX, pafsizeY, pDataSetRed->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetNIR->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, NirData[i],
		pafsizeX, pafsizeY, pDataSetNIR->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);
	pDataSetCM->RasterIO(GF_Read, offx, offy, pafsizeX, pafsizeY, cloudData[i],
		pafsizeX, pafsizeY, pDataSetCM->GetRasterBand(1)->GetRasterDataType(), 1, NULL, 0, 0, 0, 0);

	GDALClose((GDALDatasetH)pDataSetSA);
	GDALClose((GDALDatasetH)pDataSetSZ);
	GDALClose((GDALDatasetH)pDataSetVA);
	GDALClose((GDALDatasetH)pDataSetVZ);
	GDALClose((GDALDatasetH)pDataSetBlue);
	GDALClose((GDALDatasetH)pDataSetRed);
	GDALClose((GDALDatasetH)pDataSetNIR);
	GDALClose((GDALDatasetH)pdataset);
	return true;
}

bool VI::getMemory(int fileNum, int pafsizeX, int pafsizeY)
{
	BlueData = (unsigned short**)CPLMalloc(sizeof(unsigned short*)*fileNum);
	RedData = (unsigned short**)CPLMalloc(sizeof(unsigned short*)*fileNum);
	NirData = (unsigned short**)CPLMalloc(sizeof(unsigned short*)*fileNum);
	solarZenithData = (short**)CPLMalloc(sizeof(short*)*fileNum);
	viewZenithData = (short**)CPLMalloc(sizeof(short*)*fileNum);
	solarAzimuthData = (short**)CPLMalloc(sizeof(short*)*fileNum);
	viewAzimuthData = (short**)CPLMalloc(sizeof(short*)*fileNum);
	cloudData = (unsigned char**)CPLMalloc(sizeof(unsigned char*)*fileNum);
	for (int i = 0; i < fileNum; i++)
	{
		BlueData[i] = (unsigned short*)CPLMalloc(sizeof(unsigned short)*pafsizeX*pafsizeY);
		RedData[i] = (unsigned short*)CPLMalloc(sizeof(unsigned short)*pafsizeX*pafsizeY);
		NirData[i] = (unsigned short*)CPLMalloc(sizeof(unsigned short)*pafsizeX*pafsizeY);
		solarZenithData[i] = (short*)CPLMalloc(sizeof(short)*pafsizeX*pafsizeY);
		viewZenithData[i] = (short*)CPLMalloc(sizeof(short)*pafsizeX*pafsizeY);
		solarAzimuthData[i] = (short*)CPLMalloc(sizeof(short)*pafsizeX*pafsizeY);
		viewAzimuthData[i] = (short*)CPLMalloc(sizeof(short)*pafsizeX*pafsizeY);
		cloudData[i] = (unsigned char*)CPLMalloc(sizeof(unsigned char)*pafsizeX*pafsizeY);
	}
	return true;
}

bool VI::freeMemory(int fileNum)
{

	for (int i = 0; i < fileNum; i++)
	{
		CPLFree(BlueData[i]);
		CPLFree(RedData[i]);
		CPLFree(NirData[i]);
		CPLFree(solarAzimuthData[i]);
		CPLFree(solarZenithData[i]);
		CPLFree(viewAzimuthData[i]);
		CPLFree(viewZenithData[i]);
		CPLFree(cloudData[i]);
	}
	CPLFree(BlueData);
	CPLFree(RedData);
	CPLFree(NirData);
	CPLFree(solarAzimuthData);
	CPLFree(solarZenithData);
	CPLFree(viewAzimuthData);
	CPLFree(viewZenithData);
	CPLFree(cloudData);
	return true;

}
string VI::getGridNum(char * OutFileName)
{
	string namePath = (string)OutFileName;
	string name = namePath.substr(namePath.find_last_of("/") + 1);
	string s = name.substr(name.find_last_of("H"), 6);
	return s;
}

string VI::getLandCoverData(string GridNumStr, string classPath)
{
	return string();
}



vector<string> VI::getMCDPath(const char * ss)
{
	string strFileNamess(ss);
	char* strFileName = const_cast<char*> (strFileNamess.c_str());
	char* filename_str;
	filename_str = strtok(strFileName, ", ");
	vector<string> filename;
	while (filename_str != NULL)
	{
		filename.push_back(filename_str);
		filename_str = strtok(NULL, ",");
	}
	return filename;
}

float VI::mean(vector<float> data)
{
	float sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += data[i];
	}
	return sum / data.size();
}

float VI::std(vector<float> data)
{
	float avg = mean(data);
	float sum = 0.0;
	for (int i = 0; i < data.size(); i++)
	{
		sum += (data[i] - avg)*(data[i] - avg);
	}
	return sqrt(sum / data.size());
}

int VI::getMax(vector<float> data)
{
	float max = -999;
	int point = 0;
	for (int i = 0; i < data.size(); i++)
	{
		if (data[i] > max)
		{
			max = data[i];
			point = i;
		}
	}
	return point;
}

void VI::readLandCover(char * LandCoverDataPath, unsigned char * LandCoverData)
{

}

const char * VI::getFilaNameFromPath(const char * strFullPath)
{
	string file_tmp_name = strFullPath;
	string file_name = file_tmp_name.substr(file_tmp_name.find_last_of("/") + 1, file_tmp_name.find_last_of(".") - file_tmp_name.find_last_of("/") - 1);
	return file_name.c_str();
}

vector<string> VI::getMCDPath(const char * ss, const char * splitC)
{
	string strFileNamess(ss);
	char* strFileName = const_cast<char*> (strFileNamess.c_str());
	char* filename_str;
	filename_str = strtok(strFileName, splitC);
	vector<string> filename;
	while (filename_str != NULL)
	{
		filename.push_back(filename_str);
		filename_str = strtok(NULL, splitC);
	}
	return filename;
}
