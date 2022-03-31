#include <iostream>
#include <sstream>
#include "mesh.h"
#include "programlog.h"
#include "fileparser.h"

using namespace programlog;
using namespace fileparser;

namespace meshspace
{
	template <class DataType>
	bool Storage::readAny(DataType* data, int& count, string path, string paramNameInFile, string ofParamName, string fileCutStorageName)
	{
		string fileParamName = path + "inf2tr.dat";
		if (fileParser::getParam(fileParamName, paramNameInFile, count))
		{
			writeErr(PARAM_READ_ERR + " - ���������� " + ofParamName + " �� ����� " + fileParamName);
			return false;
		}

		data = new DataType[count];
		string fileStorageName = path + fileCutStorageName;
		ifstream binfile_in(fileStorageName, ios::binary);

		if (!binfile_in.is_open())
		{
			writeErr(READ_ERR + " - " + ofParamName + " - " + FILEPARSER_OPEN_FILE_ERR + " - " + fileStorageName);
			return false;
		}

		for (int i = 0; i < count; i++)
			if (!readData(binfile_in, i))
			{
				ostringstream stream_num;
				stream_num << i;
				writeErr(READ_ERR + " " + ofParamName + " - ������ � ������ " + stream_num.str() + " ���� ����������, ���� ����������� ������ � ����� " + fileStorageName);
				binfile_in.close();
				return false;
			}
		
		binfile_in.close();
		return true;
	}

	bool FinalElementMaterialStorage::readData(ifstream& binfile_in, int i)
	{
		binfile_in.read((char*)&(finalElemMaterials[i]), sizeof(int));
		if (binfile_in.gcount() != sizeof(int)) return false;
	}
	bool KnotsStorage::readData(ifstream& binfile_in, int i)
	{
		binfile_in.read((char*)&(knots[i].x), sizeof(double));
		if (binfile_in.gcount() != sizeof(double)) return false;

		binfile_in.read((char*)&(knots[i].y), sizeof(double));
		if (binfile_in.gcount() != sizeof(double)) return false;
	}

	bool FinalElementStorage::readData(ifstream& binfile_in, int i)
	{
		FinalElement endElementBuf;
		binfile_in.read((char*)&(endElementBuf.ver1), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ver2), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ver3), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ver4), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ign0), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;
		binfile_in.read((char*)&(endElementBuf.ign1), sizeof(int)); if (binfile_in.gcount() != sizeof(int)) return false;

		endElements[i] = endElementBuf;

		return true;
	}

	bool KnotsWithFirstConditionStorage::readData(ifstream& binfile_in, int i)
	{
		binfile_in.read((char*)&(IKnots[i]), sizeof(int));
		if (binfile_in.gcount() != sizeof(int)) false;

		return true;
	}

	bool MaterialStorage::readData(ifstream& in, int i)
	{
		Material materialBuf;
		in >> materialBuf.num;  if (in.failbit) return false;
		in >> materialBuf.mu;   if (in.failbit) return false;
		in >> materialBuf.toku; if (in.failbit) return false;

		
		return true;
	}

	bool FinalElementMaterialStorage::read(string path)
	{
		readAny(finalElemMaterials, ktr, path, "ktr", "������� ����������", "nvkat2d.dat");
		return true;
	}

	bool KnotsStorage::read(string path)
	{
		readAny(knots, kuzlov, path, "kuzlov", "�����", "rz.dat");
		return true;
	}

	bool FinalElementStorage::read(string path)
	{
		readAny(endElements, ktr, path, "ktr", "�������� ���������", "nvtr.dat");
		return true;
	}

	bool KnotsWithFirstConditionStorage::read(string path)
	{
		readAny(IKnots, kt1, path, "kt1", "����� � ������� �������� ���������", "l1.dat");
		return true;
	}

	bool MaterialStorage::read(string path)
	{
		readAny(materials, nMat, path, "nMat", "����������", "mat.dat");
		ifstream muFile;
		for (int i = 0; i < nMat; i++)
		{
			ostringstream stream_num;
			stream_num << materials[i].num;
			muFile.open(path + "mu.00" + stream_num.str());
			if (!muFile.is_open())
				materials[i].isMuConst = true;
			else
			{
				muFile.close();
				materials[i].isMuConst = false;
				materials[i].muSplain.init(path + "mu.00" + stream_num.str());
			}
		}
		
		return true;
	}

	Material MaterialStorage::findMaterialByNum(int num)
	{
		for (int i = 0; i < nMat; i++)
			if (materials[i].num == nMat)
				return materials[i];

		Material materialBuf;
		materialBuf.num = -1;
		return materialBuf;
	}

	bool FilePreparation::prepareMaterialFile(string path)
	{
		string numMatFile = path + "mu";
		int nMat = fileParser::getNumLinesInFile(numMatFile);
		if (nMat < 0)
			writeErr("������ ��� ���������� ������ - �� ������� ���������� ���������� ���������� ��-�� ������� � ������ - " + numMatFile);

		ofstream paramFile;
		paramFile.open(path + "nvtr.dat", ios::app);
		paramFile << " nMat= " << nMat;
		paramFile.close();

		string readMatErr = READ_ERR + " ���������� - �� ������� ���������";
		Material materialBuf;
		ifstream inMu, inToku;

		string muFileName = path + "mu";
		string tokuFileName = path + "toku";
		string materialFileName = path + "mat.dat";

		inMu.open(muFileName);
		if (!inMu.is_open())
		{
			writeErr(READ_ERR + " - mu - " + FILEPARSER_OPEN_FILE_ERR + " - " + muFileName);
			return false;
		}

		inToku.open(tokuFileName);
		if (!inMu.is_open())
		{
			writeErr(WRITE_ERR + " - " + FILEPARSER_OPEN_FILE_ERR + " - " + muFileName);
			return false;
		}

		ofstream matFile;
		matFile.open(materialFileName);


		for (int i = 0; i < nMat; i++)
		{
			inMu >> materialBuf.num;
			if (inMu.failbit)
			{
				writeMaterialErrOfReading(muFileName, "�����", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			inToku >> materialBuf.num;
			if (inToku.failbit)
			{
				writeMaterialErrOfReading(tokuFileName, "�����", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			inMu >> materialBuf.mu;
			if (inMu.failbit)
			{
				writeMaterialErrOfReading(muFileName, "mu", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			inToku >> materialBuf.toku;
			if (inToku.failbit)
			{
				writeMaterialErrOfReading(tokuFileName, "toku", i);
				closeStreams(inMu, inToku, matFile);
				return false;
			}

			matFile << materialBuf.num << " " << materialBuf.mu << " " << materialBuf.toku << endl;
			if (matFile.failbit)
			{
				writeErr(WRITE_ERR + " ��� ���������� ������ - " + FILEPARSER_OPEN_FILE_ERR + " - " + muFileName);
				closeStreams(inMu, inToku, matFile);
				return false;
			}
		}

	}

	void FilePreparation::writeMaterialErrOfReading(string fileName, string errorObj, int iErrorString)
	{
		ostringstream stream_num;
		stream_num << iErrorString;
		writeErr(READ_ERR + " ���������� ��� ���������� ������ - �� ������� ��������� " + errorObj + " ��������� � ����� " + fileName + " � ������ " + stream_num.str());
	}

	void FilePreparation::closeStreams(ifstream& in1, ifstream& in2, ofstream& out)
	{
		in1.close();
		in2.close();
		out.close();
	}

	//�������� ����������� ��������!
	Mesh::Mesh(string path)
	{
		FilePreparation::prepareMaterialFile(path);
		finalElementMaterialStorage.read(path);
		finalElementStorage.read(path);
		knotsWithFirstConditionStorage.read(path);
		materialStorage.read(path);
		knotsStorage.read(path);
	}

	Mesh::Mesh() {};
}

