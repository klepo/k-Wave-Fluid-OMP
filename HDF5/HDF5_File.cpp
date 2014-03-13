/**
 * @file        HDF5_File.cpp
 * @author      Jiri Jaros              \n
 *              CECS, ANU, Australia     \n
 *              jiri.jaros@anu.edu.au
 * 
 * @brief       The implementation file containing the HDF5 related classes
 * 
 * @version     kspaceFirstOrder3D 2.14
 * @date        27 July     2012, 14:14      (created) \n
 *              14 February 2014, 11:00 (revised
 *  
 
 * @section License
 * This file is part of the C++ extension of the k-Wave Toolbox (http://www.k-wave.org).\n
 * Copyright (C) 2012 Jiri Jaros and Bradley Treeby
 * 
 * This file is part of k-Wave. k-Wave is free software: you can redistribute it 
 * and/or modify it under the terms of the GNU Lesser General Public License as 
 * published by the Free Software Foundation, either version 3 of the License, 
 * or (at your option) any later version.
 * 
 * k-Wave is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU Lesser General Public License for more details. 
 * 
 * You should have received a copy of the GNU Lesser General Public License 
 * along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <iostream>
#include <stdexcept>
#include <unistd.h>

#include <HDF5/HDF5_File.h>

#include <Parameters/Parameters.h>
#include <Utils/ErrorMessages.h>

//----------------------------------------------------------------------------//
//-----------------------------    Constants----------------------------------//
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//

const char * THDF5_File::HDF5_MatrixDomainTypeName    = "domain_type";

const char * THDF5_File::HDF5_MatrixDataTypeName      = "data_type";

const string THDF5_File::HDF5_MatrixDomainTypeNames[] = {"real","complex"};

const string THDF5_File::HDF5_MatrixDataTypeNames[]   = {"float","long"};


const string THDF5_FileHeader::HDF5_FileTypesNames[]  = {"input","output", "checkpoint", "unknown"};

const string THDF5_FileHeader::HDF5_MajorFileVersionsNames[] = {"1"};
    
const string THDF5_FileHeader::HDF5_MinorFileVersionsNames[] = {"0","1"};

//----------------------------------------------------------------------------//
//----------------------------    THDF5_File    ------------------------------//
//------------------------------    Public     -------------------------------//
//----------------------------------------------------------------------------//


/**
 * Constructor
 */
THDF5_File::THDF5_File() :
        HDF5_FileId(H5I_BADID), FileName("")
{
    
}// end of constructor
//------------------------------------------------------------------------------




/**
 * Create the HDF5 file.
 * @param [in] FileName - File name
 * @param [in] Flags    - Flags for the HDF5 runtime
 * @throw ios:failure if error happened
 *  
 */
void THDF5_File::Create(const char * FileName, unsigned int Flags){

    // file is opened
    if (IsOpened()) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_FileCannotRecreated,FileName);        
        throw ios::failure(ErrorMessage);        
    }

    
   // Create a new file using default properties.
    this->FileName = FileName;
    
    HDF5_FileId = H5Fcreate(FileName, Flags, H5P_DEFAULT, H5P_DEFAULT);
    
    
    if (HDF5_FileId < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_FileNotCreated,FileName);        
        throw ios::failure(ErrorMessage);
    } 
    
}// end of Create
//------------------------------------------------------------------------------


/**
 * Open the HDF5 file.
 * @param [in] FileName
 * @param [in] Flags    - flags for the HDF5 runtime
 * @throw ios:failure if error happened
 *  
 */
void THDF5_File::Open(const char * FileName, unsigned int Flags){
          
    if (IsOpened())  {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_FileCannotReopen,FileName);        
        throw ios::failure(ErrorMessage);        
    };

    
    this->FileName = FileName;
    
    if (H5Fis_hdf5(FileName) == 0){
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_NotHDF5File,FileName);       
        throw ios::failure(ErrorMessage);
        
    }
    
    HDF5_FileId = H5Fopen( FileName, Flags, H5P_DEFAULT );
          
    if (HDF5_FileId < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_FileNotOpened,FileName);        
        throw ios::failure(ErrorMessage);
    }
        
    
}// end of Open
//------------------------------------------------------------------------------

/**
 * Close the HDF5 file.
 */
void THDF5_File::Close(){
         
    // Terminate access to the file.
    herr_t status = H5Fclose(HDF5_FileId);     
    if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_FileNotClosed,FileName.c_str());
        
        throw ios::failure(ErrorMessage);
    } 

    FileName    = "";
    HDF5_FileId = H5I_BADID;
    
}// end of Close
//------------------------------------------------------------------------------    

/**
 * Destructor
 */
THDF5_File::~THDF5_File(){
    
  if (IsOpened()) Close();
    
}//end of ~THDF5_File
//------------------------------------------------------------------------------


/**
 * Crate a HDF5 group.
 * @param [in] Parent  - Where to link the group at
 * @param [in] GroupName - Group name
 * @return a hanlde to the new group
 */
hid_t THDF5_File::CreateGroup(const hid_t Parent, const char * GroupName)
{
  hid_t HDF5_group_id = H5Gcreate(Parent, GroupName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  //if error
  if (HDF5_group_id == H5I_INVALID_HID)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage,HDF5_ERR_FMT_GroupNotCreated,GroupName, FileName.c_str());        
    throw ios::failure(ErrorMessage);    
  }
  
  return HDF5_group_id;
};// end of CreateGroup
//------------------------------------------------------------------------------


/**
 * Open a HDF5 group
 * @param [in] ParentGroupOrFile
 * @param [in] GroupName
 * @return 
 */
hid_t THDF5_File::OpenGroup(const hid_t Parent, const char * GroupName)
{
  hid_t 
  HDF5_group_id = H5Gopen(Parent, GroupName, H5P_DEFAULT);
  
  //if error
  if (HDF5_group_id == H5I_INVALID_HID)
  {
    char ErrorMessage[256];
    sprintf(ErrorMessage,HDF5_ERR_FMT_GroupNotOpened,GroupName, FileName.c_str());        
    throw ios::failure(ErrorMessage);    
  }
  
  return HDF5_group_id;  
}// end of OpenGroup
//------------------------------------------------------------------------------


/**
 * Close a group
 * @param[in] Group
 */
void THDF5_File::CloseGroup(const hid_t HDF5_group_id)
{  
  H5Gclose(HDF5_group_id);  
}// end of CloseGroup
//------------------------------------------------------------------------------

/**
 * Open the dataset.
 * @throw ios::failure
 * @param [in]  DatasetName
 * @return      Dataset id
 * 
 */
hid_t THDF5_File::OpenDataset  (const char * DatasetName){
  
      // Open dataset
    hid_t HDF5_dataset_id = H5Dopen(HDF5_FileId, DatasetName, H5P_DEFAULT);
    
    if (HDF5_dataset_id == H5I_INVALID_HID ){
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_DatasetNotOpened,FileName.c_str(), DatasetName);        
        throw ios::failure(ErrorMessage);    
    }
    
    return HDF5_dataset_id;     
    
}// end of OpenDataset
//------------------------------------------------------------------------------


/**
 * Create the dataset. 
 * @throw ios::failure
 * 
 * @param [in]  DatasetName - Dataset name         
 * @param [in]  DimensionSizes - Dataset dimension sizes
 * @param [in]  ChunkSizes  - 3D Chunk size
 * @param [in]  CompressionLevel - Compression level
 * @return      Dataset_id
 */ 
hid_t THDF5_File::CreateFloatDataset(const char * DatasetName, const TDimensionSizes & DimensionSizes, const TDimensionSizes & ChunkSizes, const int CompressionLevel){
 
    const int RANK = 3;
    hsize_t    Dims [RANK] = {DimensionSizes.Z, DimensionSizes.Y, DimensionSizes.X};
    hsize_t    Chunk[RANK] = {ChunkSizes.Z, ChunkSizes.Y, ChunkSizes.X};
    hid_t      Property_list;
    herr_t     Status;
    
    hid_t dataspace_id = H5Screate_simple (RANK, Dims, NULL); 
    
    
    // set chunk size
    Property_list  = H5Pcreate (H5P_DATASET_CREATE);
    
    Status = H5Pset_chunk(Property_list, RANK, Chunk);    
    if (Status < 0 ){
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_DatasetNotOpened,FileName.c_str(), DatasetName);        
        throw ios::failure(ErrorMessage);    
       
    }
    
    // set compression level
    Status = H5Pset_deflate (Property_list, CompressionLevel); 
    if (Status < 0 ){
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotSetCompression,FileName.c_str(), DatasetName, CompressionLevel);        
        throw ios::failure(ErrorMessage);    
       
    }
    
    // create dataset
    hid_t HDF5_dataset_id = H5Dcreate (HDF5_FileId, DatasetName, H5T_NATIVE_FLOAT, dataspace_id,
                            H5P_DEFAULT, Property_list, H5P_DEFAULT);
    
    if (HDF5_dataset_id == H5I_INVALID_HID ){
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_DatasetNotOpened,FileName.c_str(), DatasetName);        
        throw ios::failure(ErrorMessage);    
       
    }
    
    H5Pclose (Property_list);
    
    return HDF5_dataset_id; 
    
}// end of CreateDataset    
//------------------------------------------------------------------------------




/**
 * Close dataset.
 * @param [in] HDF5_Dataset_id
 *  
 */
void  THDF5_File::CloseDataset(const hid_t HDF5_Dataset_id)
{         
  H5Dclose (HDF5_Dataset_id);    
}// end of CloseDataset
//------------------------------------------------------------------------------




/**
 * Write data into a dataset.
 * 
 * @param [in] DatasetName
 * @param [in] DimensionSizes
 * @param [in] Data
 * @throw ios::failure
 */
void THDF5_File::WriteCompleteDataset(const char * DatasetName, const TDimensionSizes & DimensionSizes, const float * Data){
  
  const int rank = 3;
  const hsize_t dims[]={DimensionSizes.Z, DimensionSizes.Y,DimensionSizes.X};
  
  // write to dataset
  herr_t  status = H5LTmake_dataset_float(HDF5_FileId,DatasetName,rank,dims, Data);
  
   if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotWriteTo,DatasetName);
        
        throw ios::failure(ErrorMessage);
    } 
  
}// end of WriteDataset (float)
//------------------------------------------------------------------------------


/**
 * Write data into a dataset.
 * 
 * @param [in] DatasetName
 * @param [in] DimensionSizes
 * @param [in] Data
 * @throw ios::failure
 */
void THDF5_File::WriteCompleteDataset (const char * DatasetName, const TDimensionSizes & DimensionSizes, const long  * Data){
  const int     rank = 3;
  const hsize_t dims[]={DimensionSizes.Z, DimensionSizes.Y,DimensionSizes.X};
    
  // write to dataset
  herr_t  status = H5LTmake_dataset(HDF5_FileId,DatasetName,rank,dims, H5T_STD_U64LE, Data);
  
   if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotWriteTo,DatasetName);
        
        throw ios::failure(ErrorMessage);
    } 
    
}// end of WriteDataset (long)
//------------------------------------------------------------------------------



/**
 * Write a hyperslab into the dataset.
 * 
 * @param [in] HDF5_Dataset_id - Dataset id
 * @param [in] Position - Position in the dataset
 * @param [in] Size - Size of the hyperslab
 * @param [in] Data
 * @throw ios::failure
 * 
 */
void THDF5_File::WriteHyperSlab(const hid_t HDF5_Dataset_id, const TDimensionSizes & Position , const TDimensionSizes & Size, const float * Data){
    
    
    
     // Select hyperslab
    const int MatrixRank = 3;    
    hsize_t ElementCount[MatrixRank] = {Size.Z, Size.Y, Size.X};    
    hsize_t Offset[MatrixRank] = {Position.Z,Position.Y,Position.X};
    
    herr_t status;
    hid_t  HDF5_Filespace,HDF5_Memspace;
    
    // Select hyperslab in the file.     
    HDF5_Filespace = H5Dget_space(HDF5_Dataset_id);  
    
    
    status = H5Sselect_hyperslab(HDF5_Filespace, H5S_SELECT_SET, Offset, 0, ElementCount, NULL);
    if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotWriteTo,"");        
        throw ios::failure(ErrorMessage);
     } 
    
    
    // assign memspace
    HDF5_Memspace = H5Screate_simple(MatrixRank, ElementCount, NULL);

    
     status = H5Dwrite(HDF5_Dataset_id, H5T_NATIVE_FLOAT, HDF5_Memspace, HDF5_Filespace,  H5P_DEFAULT, Data);
     if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotWriteTo,"");
        
        throw ios::failure(ErrorMessage);
     } 
    
    H5Sclose(HDF5_Memspace);
    H5Sclose(HDF5_Filespace);  
    
    
    
}// end of WriteHyperSlab
//------------------------------------------------------------------------------


/**
 * Write hyperslab.
 * 
 * @param [in] HDF5_Dataset_id - Dataset id
 * @param [in] Position - Position in the dataset
 * @param [in] Size - Size of the hyperslab
 * @param [in] Data
 * @throw ios::failure
 */
void THDF5_File::WriteHyperSlab(const hid_t HDF5_Dataset_id, const TDimensionSizes & Position , const TDimensionSizes & Size, const long * Data){
    
    
    
     // Select hyperslab
    const int MatrixRank = 3;    
    hsize_t ElementCount[MatrixRank] = {Size.Z, Size.Y, Size.X};    
    hsize_t Offset[MatrixRank] = {Position.Z,Position.Y,Position.X};
    
    herr_t status;
    hid_t  HDF5_Filespace,HDF5_Memspace;
    
    // Select hyperslab in the file.     
    HDF5_Filespace = H5Dget_space(HDF5_Dataset_id);  
    
    
    status = H5Sselect_hyperslab(HDF5_Filespace, H5S_SELECT_SET, Offset, 0, ElementCount, NULL);
    if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotWriteTo,"");        
        throw ios::failure(ErrorMessage);
     } 
    
    
    // assign memspace
    HDF5_Memspace = H5Screate_simple(MatrixRank, ElementCount, NULL);

    
     status = H5Dwrite(HDF5_Dataset_id, H5T_STD_U64LE, HDF5_Memspace, HDF5_Filespace,  H5P_DEFAULT, Data);
     if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotWriteTo,"");
        
        throw ios::failure(ErrorMessage);
     } 
    
    H5Sclose(HDF5_Memspace);
    H5Sclose(HDF5_Filespace);  
    
    
    
}// end of WriteHyperSlab
//------------------------------------------------------------------------------



/**
 * Write a scalar value into the file.
 * @param [in] DatasetName
 * @param [in] Value
 * 
 */
void THDF5_File::WriteScalarValue(const char * DatasetName, const float Value){
   
    WriteCompleteDataset (DatasetName, TDimensionSizes(1,1,1), &Value);
    WriteMatrixDataType  (DatasetName, hdf5_mdt_float);        
    WriteMatrixDomainType(DatasetName, hdf5_mdt_real);
    
    
}// end of WriteScalarValue (float)
//------------------------------------------------------------------------------

/**
 * Write a scalar value.
 * @param [in] DatasetName
 * @param [in] Value
 * 
 */
void THDF5_File::WriteScalarValue(const char * DatasetName, const long Value){
    
    WriteCompleteDataset (DatasetName, TDimensionSizes(1,1,1), &Value);
    WriteMatrixDataType  (DatasetName, hdf5_mdt_long);        
    WriteMatrixDomainType(DatasetName, hdf5_mdt_real);
    
    
    
}// end of WriteScalarValue
//------------------------------------------------------------------------------


/**
 * Read data from the dataset.
 * 
 * @param [in]  DatasetName
 * @param [in]  DimensionSizes
 * @param [out] Data
 * @throw ios::failure
 */
void THDF5_File::ReadCompleteDataset (const char * DatasetName, const TDimensionSizes & DimensionSizes, float * Data){
   
    if (GetDatasetDimensionSizes(DatasetName).GetElementCount() != 
                               DimensionSizes.GetElementCount()) {    
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_WrongDimensionSizes,DatasetName);        
        throw ios::failure(ErrorMessage);        
    }

     /* read dataset */
    herr_t status = H5LTread_dataset_float(HDF5_FileId,DatasetName,Data);
    if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotReadFrom,DatasetName);
        
        throw ios::failure(ErrorMessage);
    } 
 
}// end of ReadDataset (float)
//------------------------------------------------------------------------------

/**
 * Write data from the dataset.
 * 
 * @param [in]  DatasetName
 * @param [in]  DimensionSizes
 * @param [out] Data
 * @throw ios::failure
 */
void THDF5_File::ReadCompleteDataset  (const char * DatasetName, const TDimensionSizes & DimensionSizes,long  * Data){

    if (GetDatasetDimensionSizes(DatasetName).GetElementCount() != 
                               DimensionSizes.GetElementCount()) {    
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_WrongDimensionSizes,DatasetName);        
        throw ios::failure(ErrorMessage);       
    }           
            
     /* read dataset */
    herr_t status = H5LTread_dataset(HDF5_FileId,DatasetName,H5T_STD_U64LE, Data);
    if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotReadFrom,DatasetName);
        
        throw ios::failure(ErrorMessage);
    }   
}// end of ReadDataset(long)
//------------------------------------------------------------------------------



/**
 * Get dimension sizes of the dataset.
 * @param [in] DatasetName
 * @return     DimensionSizes
 * @throw ios::failure
 */
 TDimensionSizes THDF5_File::GetDatasetDimensionSizes(const char * DatasetName){
  
     hsize_t dims[3] = {0,0,0};
            
     herr_t status = H5LTget_dataset_info(HDF5_FileId,DatasetName,dims,NULL,NULL);    
     if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotReadFrom,DatasetName);
        
        throw ios::failure(ErrorMessage);
     } 
     
     return TDimensionSizes(dims[2],dims[1],dims[0]);
          
 }// end of CheckDatasetDimensions
 //------------------------------------------------------------------------------
 

/**
 * Get number of elements of the dataset.
 * @param [in] DatasetName
 * @return     Number of elements
 * @throw ios::failure
 */
size_t THDF5_File::GetDatasetElementCount(const char * DatasetName){

    hsize_t dims[3] = {0,0,0};
            
     herr_t status = H5LTget_dataset_info(HDF5_FileId,DatasetName,dims,NULL,NULL);    
     if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotReadFrom,DatasetName);
        
        throw ios::failure(ErrorMessage);
     } 
     
     return dims[0] * dims[1] * dims[2];
 }// end of GetDatasetElementCount
 //-----------------------------------------------------------------------------
 

 
 
/**
 * Write matrix data type into the dataset.
 * 
 * @param [in] DatasetName
 * @param [in] MatrixDataType
 */
void THDF5_File::WriteMatrixDataType(const char * DatasetName, const THDF5_MatrixDataType & MatrixDataType){
    
   WriteStringAttribute(DatasetName, HDF5_MatrixDataTypeName, HDF5_MatrixDataTypeNames[MatrixDataType]);
    
}// end of WriteMatrixDataType
//------------------------------------------------------------------------------


/**
 * Write matrix data type into the dataset.
 * 
 * @param [in] DatasetName
 * @param [in] MatrixDomainType
 */
void THDF5_File::WriteMatrixDomainType(const char * DatasetName, const THDF5_MatrixDomainType & MatrixDomainType){
       
    WriteStringAttribute(DatasetName, HDF5_MatrixDomainTypeName, HDF5_MatrixDomainTypeNames [MatrixDomainType]);
        
}// end of WriteMatrixDomainType
//------------------------------------------------------------------------------    


/**
 * Read matrix data type.
 * 
 * @param [in]  DatasetName
 * @return      MatrixDataType
 * @throw ios::failure
 */    
THDF5_File::THDF5_MatrixDataType THDF5_File::ReadMatrixDataType(const char * DatasetName){
        
    string ParamValue;
    
    ParamValue = ReadStringAttribute(DatasetName, HDF5_MatrixDataTypeName);
    
    if (ParamValue == HDF5_MatrixDataTypeNames[0]) return (THDF5_MatrixDataType) 0; 
    if (ParamValue == HDF5_MatrixDataTypeNames[1]) return (THDF5_MatrixDataType) 1;
    
    char ErrorMessage[256];
    sprintf(ErrorMessage,HDF5_ERR_FMT_BadAttributeValue,DatasetName, HDF5_MatrixDataTypeName,ParamValue.c_str());        
    throw ios::failure(ErrorMessage);
    
    return (THDF5_MatrixDataType) 0; 
    
}// end of ReadMatrixDataType
//------------------------------------------------------------------------------


/**
 * Read matrix domain type.
 * 
 * @param [in]  DatasetName
 * @return      DomainType
 * @throw ios::failure
 */    
THDF5_File::THDF5_MatrixDomainType THDF5_File::ReadMatrixDomainType(const char * DatasetName){
        
    string ParamValue;
         
    ParamValue = ReadStringAttribute(DatasetName, HDF5_MatrixDomainTypeName);    
    
    if (ParamValue == HDF5_MatrixDomainTypeNames[0]) return (THDF5_MatrixDomainType) 0; 
    if (ParamValue == HDF5_MatrixDomainTypeNames[1]) return (THDF5_MatrixDomainType) 1;
    
    
    char ErrorMessage[256];
    sprintf(ErrorMessage,HDF5_ERR_FMT_BadAttributeValue,DatasetName, HDF5_MatrixDomainTypeName,ParamValue.c_str());        
    throw ios::failure(ErrorMessage);
   
    
    return (THDF5_MatrixDomainType) 0; 
        
}// end of ReadMatrixDomainType
//------------------------------------------------------------------------------    
    

/**
 * Write integer attribute.
 * 
 * @param [in] DatasetName
 * @param [in] AttributeName
 * @param [in] Value
 * @throw ios::failure
 */
inline void THDF5_File::WriteStringAttribute(const char * DatasetName, const char * AttributeName, const string & Value){
    
    herr_t status = H5LTset_attribute_string(HDF5_FileId, DatasetName, AttributeName, Value.c_str());
    
    if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotWriteToAttribute,AttributeName, DatasetName);
        
        throw ios::failure(ErrorMessage);
     } 
    
}// end of WriteIntAttribute
//------------------------------------------------------------------------------



/**
 * Read integer attribute.
 * 
 * @param  [in] DatasetName
 * @param  [in] AttributeName
 * @return  Attribute value
 * @throw ios::failure
 */
inline  string THDF5_File::ReadStringAttribute(const char * DatasetName, const char * AttributeName){
        
    char Value[256] = "";        
    herr_t status = H5LTget_attribute_string (HDF5_FileId, DatasetName, AttributeName, Value);            
    
    if (status < 0) {
        char ErrorMessage[256];
        sprintf(ErrorMessage,HDF5_ERR_FMT_CouldNotReadFromAttribute,AttributeName, DatasetName);
        
        throw ios::failure(ErrorMessage);
    }
    
    return Value;
        
}// end of ReadIntAttribute
//------------------------------------------------------------------------------


//----------------------------------------------------------------------------//
//----------------------------    THDF5_File    ------------------------------//
//----------------------------    Protected     ------------------------------//
//----------------------------------------------------------------------------//








//----------------------------------------------------------------------------//
//------------------------    THDF5_File_Header ------------------------------//
//-----------------------------    Public     --------------------------------//
//----------------------------------------------------------------------------//

/**
 * Constructor.
 */
THDF5_FileHeader::THDF5_FileHeader(){
   
    HDF5_FileHeaderValues.clear();
    PopulateHeaderFileMap();
    
}// end of constructor
//------------------------------------------------------------------------------


/**
 * Copy constructor.
 * @param [in] other
 */
THDF5_FileHeader::THDF5_FileHeader(const THDF5_FileHeader & other){
    HDF5_FileHeaderValues = other.HDF5_FileHeaderValues;
    HDF5_FileHeaderNames  = other.HDF5_FileHeaderNames;
    
}// end of copy constructor
//------------------------------------------------------------------------------


/**
 * Destructor.
 * 
 */
THDF5_FileHeader::~THDF5_FileHeader(){
    HDF5_FileHeaderValues.clear();
    HDF5_FileHeaderNames.clear();
    
}// end of destructor
//------------------------------------------------------------------------------



/**
 * Read header from the input file.
 * @param [in] InputFile - Input file to read from
 */
void THDF5_FileHeader::ReadHeaderFromInputFile(THDF5_File & InputFile){
    
    // read file type
    HDF5_FileHeaderValues[hdf5_fhi_file_type] = InputFile.ReadStringAttribute("/",HDF5_FileHeaderNames[hdf5_fhi_file_type].c_str());
    
    if (GetFileType() == hdf5_ft_input) {                
        HDF5_FileHeaderValues[hdf5_fhi_created_by]       = InputFile.ReadStringAttribute("/",HDF5_FileHeaderNames[hdf5_fhi_created_by].c_str());
        HDF5_FileHeaderValues[hdf5_fhi_creation_date]    = InputFile.ReadStringAttribute("/",HDF5_FileHeaderNames[hdf5_fhi_creation_date].c_str());
        HDF5_FileHeaderValues[hdf5_fhi_file_description] = InputFile.ReadStringAttribute("/",HDF5_FileHeaderNames[hdf5_fhi_file_description].c_str());
        HDF5_FileHeaderValues[hdf5_fhi_major_version]    = InputFile.ReadStringAttribute("/",HDF5_FileHeaderNames[hdf5_fhi_major_version].c_str());
        HDF5_FileHeaderValues[hdf5_fhi_minor_version]    = InputFile.ReadStringAttribute("/",HDF5_FileHeaderNames[hdf5_fhi_minor_version].c_str());        
    }
    
    
}// end of ReadHeaderFromInputFile
//------------------------------------------------------------------------------

/**
 * Write header into the output file.
 * @param [in] OutputFile
 */
void THDF5_FileHeader::WriteHeaderToOutputFile (THDF5_File & OutputFile){
    
    for (map<THDF5_FileHeaderItems, string>::iterator it = HDF5_FileHeaderNames.begin(); 
         it != HDF5_FileHeaderNames.end(); it++) {
        
        OutputFile.WriteStringAttribute("/",it->second.c_str(), HDF5_FileHeaderValues[it->first]);
    }
    
}// end of WriteHeaderToOutputFile
//------------------------------------------------------------------------------


/**
 * Get File type. 
 * @return FileType
 */
THDF5_FileHeader::THDF5_FileType  THDF5_FileHeader::GetFileType(){
    
    for (int i = hdf5_ft_input; i < hdf5_ft_unknown ; i++){
      if (HDF5_FileHeaderValues[hdf5_fhi_file_type] == HDF5_FileTypesNames[static_cast<THDF5_FileType >(i)]) return static_cast<THDF5_FileType >(i);
    }
    
    return THDF5_FileHeader::hdf5_ft_unknown;
    
}// end of GetFileType
//------------------------------------------------------------------------------
    
/**
 * Set File type.
 * @param FileType
 */
void THDF5_FileHeader::SetFileType(const THDF5_FileHeader::THDF5_FileType FileType){
    
    HDF5_FileHeaderValues[hdf5_fhi_file_type] = HDF5_FileTypesNames[FileType];
    
}// end of SetFileType
//------------------------------------------------------------------------------    


/**
 * Get File Version an enum
 * @return file version as an enum
 */
THDF5_FileHeader::THDF5_FileVersion THDF5_FileHeader::GetFileVersion(){
            
    if ((HDF5_FileHeaderValues[hdf5_fhi_major_version] == HDF5_MajorFileVersionsNames[0]) &&
        (HDF5_FileHeaderValues[hdf5_fhi_minor_version] == HDF5_MinorFileVersionsNames[0]))         
       return hdf5_fv_10;
    
    if ((HDF5_FileHeaderValues[hdf5_fhi_major_version] == HDF5_MajorFileVersionsNames[0]) &&
        (HDF5_FileHeaderValues[hdf5_fhi_minor_version] == HDF5_MinorFileVersionsNames[1]))         
       return hdf5_fv_11;
        
    return hdf5_fv_unknown;
            
}// end of GetFileVersion
//------------------------------------------------------------------------------    

/**
 * Set actual date and time. 
 * 
 */
void THDF5_FileHeader::SetActualCreationTime(){
        
    struct tm *current;
    time_t now;    
    time(&now);
    current = localtime(&now);     
    
    char DateString[20];     
    
    sprintf(DateString, "%02i/%02i/%02i, %02i:%02i:%02i", 
                        current->tm_mday, current->tm_mon+1, current->tm_year-100,
                        current->tm_hour, current->tm_min, current->tm_sec);
    
    HDF5_FileHeaderValues[hdf5_fhi_creation_date] = DateString;
    
}// end of SetCreationTime
//------------------------------------------------------------------------------

/**
 * Set Host name.
 * 
 */
void THDF5_FileHeader::SetHostName(){
    
    char   HostName[256];     
    gethostname(HostName, 256);
    
    HDF5_FileHeaderValues[hdf5_fhi_host_name] = HostName;
    

}// end of SetHostName
//------------------------------------------------------------------------------


/**
 * Set memory consumption.
 * @param [in] TotalMemory
 */
void THDF5_FileHeader::SetMemoryConsumption(size_t TotalMemory){
    
    char Text[20] = "";
    sprintf(Text, "%ld MB",TotalMemory);
       
   
    HDF5_FileHeaderValues[hdf5_fhi_total_memory_consumption]     = Text;
    
    sprintf(Text, "%ld MB",TotalMemory / TParameters::GetInstance()->GetNumberOfThreads());
    HDF5_FileHeaderValues[hdf5_fhi_peak_core_memory_consumption] = Text;
    
}// end of SetMemoryConsumption
//------------------------------------------------------------------------------


/** 
 * Set execution times in file header.
 * @param [in] TotalTime  
 * @param [in] LoadTime
 * @param [in] PreProcessingTime 
 * @param [in] SimulationTime
 * @param [in] PostprocessingTime
 */
void THDF5_FileHeader::SetExecutionTimes(const double TotalTime, const double LoadTime,
                           const double PreProcessingTime, const double SimulationTime,
                           const double PostprocessingTime){
    
    char Text [30] = "";
    
    sprintf(Text,"%8.2fs", TotalTime);
    HDF5_FileHeaderValues[hdf5_fhi_total_execution_time] = Text;        
    
    sprintf(Text,"%8.2fs", LoadTime);
    HDF5_FileHeaderValues[hdf5_fhi_data_load_time] = Text;       
        
    sprintf(Text,"%8.2fs", PreProcessingTime);
    HDF5_FileHeaderValues[hdf5_fhi_preprocessing_time] = Text;       
        
    
    sprintf(Text,"%8.2fs", SimulationTime);
    HDF5_FileHeaderValues[hdf5_fhi_simulation_time] = Text;       
    
    sprintf(Text,"%8.2fs", PostprocessingTime);
    HDF5_FileHeaderValues[hdf5_fhi_postprocessing_time] = Text;       
                             
                          
}// end of SetExecutionTimes
//------------------------------------------------------------------------------



/**
 * Set Number of cores.
 * 
 */
void THDF5_FileHeader::SetNumberOfCores(){
    char Text[12] = "";
    sprintf(Text, "%d",TParameters::GetInstance()->GetNumberOfThreads());
    
    HDF5_FileHeaderValues[hdf5_fhi_number_of_cores] = Text;
    
}// end of SetNumberOfCores
//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//------------------------    THDF5_File_Header ------------------------------//
//---------------------------    Protected     --------------------------------//
//----------------------------------------------------------------------------//






/**
 * Create map with names for the header. 
 * 
 */
void THDF5_FileHeader::PopulateHeaderFileMap(){
    
    HDF5_FileHeaderNames.clear();
    
    HDF5_FileHeaderNames[hdf5_fhi_created_by]                   = "created_by";
    HDF5_FileHeaderNames[hdf5_fhi_creation_date]                = "creation_date";
    HDF5_FileHeaderNames[hdf5_fhi_file_description]             = "file_description";
    HDF5_FileHeaderNames[hdf5_fhi_major_version]                = "major_version";
    HDF5_FileHeaderNames[hdf5_fhi_minor_version]                = "minor_version";
    HDF5_FileHeaderNames[hdf5_fhi_file_type]                    = "file_type";
    
    HDF5_FileHeaderNames[hdf5_fhi_host_name]                    = "host_names";
    HDF5_FileHeaderNames[hdf5_fhi_number_of_cores]              = "number_of_cpu_cores" ;
    HDF5_FileHeaderNames[hdf5_fhi_total_memory_consumption]     = "total_memory_in_use";
    HDF5_FileHeaderNames[hdf5_fhi_peak_core_memory_consumption] = "peak_core_memory_in_use";
    
    HDF5_FileHeaderNames[hdf5_fhi_total_execution_time]         = "total_execution_time";
    HDF5_FileHeaderNames[hdf5_fhi_data_load_time]               = "data_loading_phase_execution_time";
    HDF5_FileHeaderNames[hdf5_fhi_preprocessing_time]           = "pre-processing_phase_execution_time";    
    HDF5_FileHeaderNames[hdf5_fhi_simulation_time]              = "simulation_phase_execution_time";
    HDF5_FileHeaderNames[hdf5_fhi_postprocessing_time]          = "post-processing_phase_execution_time";
                                
}// end of PopulateHeaderFileMap
//------------------------------------------------------------------------------