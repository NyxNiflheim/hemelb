// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

// LocalPropertyHdf5Output.cc
// single h5 file + Group + compound + xmf

#ifdef USE_HDF5

#include "extraction/LocalPropertyHdf5Output.h"
#include "net/IOCommunicator.h"
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <filesystem>

// HDF5错误检查宏
#define H5_CHECK(err) (h5_check((err), __FILE__, __LINE__))
inline herr_t h5_check(herr_t err, const char* file, int line) {
    if (err < 0) {
        fprintf(stderr, "HDF5 Error at %s:%d\n", file, line);
        H5Eprint(H5E_DEFAULT, stderr);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    return err;
}

namespace hemelb::extraction
{
  LocalPropertyHdf5Output::LocalPropertyHdf5Output(IterableDataSource& dataSource,
                                                     const PropertyOutputFile& outputSpec,
                                                     const net::IOCommunicator& ioComms)
    : LocalPropertyOutput(dataSource, outputSpec, comms), mpi_comm(comms) {}

  LocalPropertyHdf5Output::~LocalPropertyHdf5Output()
  {
    if (file_id >= 0) {
      WriteXDMFFile();
      if(compound_type_id >= 0) H5Tclose(compound_type_id);
      H5Fclose(file_id);
    }
  }

  void LocalPropertyHdf5Output::CreateCompoundType() {
    // Pass 1: 计算内存布局
    size_t current_offset = 0;
    field_offsets.clear();
    
    field_offsets.push_back(current_offset); current_offset += sizeof(uint32_t); // i
    field_offsets.push_back(current_offset); current_offset += sizeof(uint32_t); // j
    field_offsets.push_back(current_offset); current_offset += sizeof(uint32_t); // k
    
    for (const auto& field_spec : outputSpec.fields) {
      field_offsets.push_back(current_offset);
      // 根据typecode来累加大小
      current_offset += GetFieldLength(field_spec.src) * std::visit([](auto t){ return sizeof(decltype(t)); }, field_spec.typecode);
    }
    compound_type_size = current_offset;

    // Pass 2: 创建HDF5复合类型
    compound_type_id = H5Tcreate(H5T_COMPOUND, compound_type_size); 

    H5Tinsert(compound_type_id, "i", field_offsets[0], H5T_NATIVE_UINT32);
    H5Tinsert(compound_type_id, "j", field_offsets[1], H5T_NATIVE_UINT32);
    H5Tinsert(compound_type_id, "k", field_offsets[2], H5T_NATIVE_UINT32);
    
    for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
      const auto& field_spec = outputSpec.fields[i];
      unsigned field_len = GetFieldLength(field_spec.src);
      size_t member_offset = field_offsets[i + 3];

      // 访问variant，为每种类型选择正确的HDF5类型
      hid_t hdf_type = std::visit([](auto t) -> hid_t {
          using T = decltype(t);
          if constexpr (std::is_same_v<T, int>) return H5T_NATIVE_INT;
          if constexpr (std::is_same_v<T, float>) return H5T_NATIVE_FLOAT;
          if constexpr (std::is_same_v<T, double>) return H5T_NATIVE_DOUBLE;
          return H5T_NATIVE_VOID; // Should not happen
      }, field_spec.typecode);

      if (field_len == 1) {
        H5Tinsert(compound_type_id, field_spec.name.c_str(), member_offset, hdf_type);
      } else {
        hsize_t dims[1] = {field_len};
        hid_t array_tid = H5Tarray_create(hdf_type, 1, dims);
        H5Tinsert(compound_type_id, field_spec.name.c_str(), member_offset, array_tid);
        H5Tclose(array_tid);
      }
    }
  }

  void LocalPropertyHdf5Output::WriteXDMFFile() {
    if (comms.Rank() != 0) return;
    
    std::string h5_filename = outputSpec.filename.string();
    std::string xmf_filename = h5_filename;
    size_t pos = xmf_filename.rfind(".h5");
    if (pos != std::string::npos) xmf_filename.replace(pos, 3, ".xmf");
    else xmf_filename += ".xmf";

    std::filesystem::path h5_path(h5_filename);
    std::string h5_basename = h5_path.filename().string();
    
    std::ofstream xmf_file(xmf_filename);
    xmf_file << "<?xml version=\"1.0\" ?>\n";
    xmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    xmf_file << "<Xdmf Version=\"3.0\">\n";
    xmf_file << "  <Domain>\n";
    xmf_file << "    <Grid Name=\"HemeLB Collection\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

    for (unsigned long ts : written_timesteps) {
        xmf_file << "      <Grid Name=\"step_" << ts << "\" GridType=\"Uniform\">\n";
        xmf_file << "        <Time Value=\"" << ts << "\"/>\n";
        xmf_file << "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << global_site_count << "\"/>\n";
        
        std::string hdf5_path = h5_basename + ":/step_" + std::to_string(ts) + "/results";

        xmf_file << "        <Geometry GeometryType=\"XYZ\">\n";
        xmf_file << "          <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << global_site_count << " 3\" Type=\"HyperSlab\">\n";
        xmf_file << "            <DataItem Dimensions=\"3 2\" Format=\"XML\">0 0 1 " << global_site_count << " 1 3</DataItem>\n";
        xmf_file << "            <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " " << compound_type_size << "\" DataType=\"Char\">" << hdf5_path << "</DataItem>\n";
        xmf_file << "          </DataItem>\n";
        xmf_file << "        </Geometry>\n";

        for (const auto& field : outputSpec.fields) {
            unsigned field_len = GetFieldLength(field.src);
            std::string type = (field_len > 1) ? "Vector" : "Scalar";
            xmf_file << "        <Attribute Name=\"" << field.name << "\" AttributeType=\"" << type << "\" Center=\"Node\">\n";
            xmf_file << "          <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << global_site_count << " " << field_len << "\" Type=\"HyperSlab\">\n";
            xmf_file << "            <DataItem Dimensions=\"3 2\" Format=\"XML\">0 " << (&field - &outputSpec.fields[0]) + 3 << " 1 " << global_site_count << " 1 " << field_len << "</DataItem>\n";
            xmf_file << "            <DataItem Format=\"HDF\" Dimensions=\"" << global_site_count << " " << compound_type_size << "\" DataType=\"Char\">" << hdf5_path << "</DataItem>\n";
            xmf_file << "          </DataItem>\n";
            xmf_file << "        </Attribute>\n";
        }
        xmf_file << "      </Grid>\n";
    }

    xmf_file << "    </Grid>\n";
    xmf_file << "  </Domain>\n";
    xmf_file << "</Xdmf>\n";
    xmf_file.close();
  }

  void LocalPropertyHdf5Output::Write(unsigned long timestepNumber, unsigned long totalSteps) {
    if (!ShouldWrite(timestepNumber)) return;
    
    if (file_id < 0) { // first write, create a new file
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5_CHECK(H5Pset_fapl_mpio(fapl, mpi_comm, MPI_INFO_NULL));
        file_id = H5Fcreate(outputSpec.filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5_CHECK(H5Pclose(fapl));
        CreateCompoundType();
    } 
    // 如果不是第一次写，需要以读写模式重新打开文件
    if (!written_timesteps.empty()) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5_CHECK(H5Pset_fapl_mpio(fapl, mpi_comm, MPI_INFO_NULL));
        file_id = H5Fopen(outputSpec.filename.c_str(), H5F_ACC_RDWR, fapl);
        H5_CHECK(H5Pclose(fapl));
    }
    
    written_timesteps.push_back(timestepNumber);

    std::string group_name = "step_" + std::to_string(timestepNumber);
    hid_t group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHECK(group_id);
    
    hsize_t global_dims[1] = {global_site_count};
    hid_t filespace_id = H5Screate_simple(1, global_dims, NULL);
    hid_t dataset_id = H5Dcreate2(group_id, "results", compound_type_id, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHECK(dataset_id);
    
    std::vector<char> buffer(local_site_count * compound_type_size);
    char* buffer_ptr = buffer.data();

    dataSource.Reset();
    while(dataSource.ReadNext()) {
      if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition())) {
        const auto& pos = dataSource.GetPosition();
        *(reinterpret_cast<uint32_t*>(buffer_ptr + field_offsets[0])) = pos.x();
        *(reinterpret_cast<uint32_t*>(buffer_ptr + field_offsets[1])) = pos.y();
        *(reinterpret_cast<uint32_t*>(buffer_ptr + field_offsets[2])) = pos.z();
        
        for (size_t i = 0; i < outputSpec.fields.size(); ++i) {
          const auto& field_spec = outputSpec.fields[i];
          char* field_ptr = buffer_ptr + field_offsets[i + 3];
          std::visit([&](auto&& src_type) {
            using T = decltype(type_tag);
            T* dest = reinterpret_cast<T*>(field_ptr);

             if constexpr (std::is_same_v<T, double>) {
                 if constexpr (std::is_same_v<std::decay_t<decltype(field_spec.src)>, source::Velocity>) {
                    const auto& vel = dataSource.GetVelocity();
                    dest[0] = std::isinf(vel.x()) ? 0.0 : vel.x();
                    dest[1] = std::isinf(vel.y()) ? 0.0 : vel.y();
                    dest[2] = std::isinf(vel.z()) ? 0.0 : vel.z();
                 }
             }else if constexpr (std::is_same_v<T, int>) {
                 if constexpr (std::is_same_v<std::decay_t<decltype(field_spec.src)>, source::MpiRank>) {
                    dest[0] = comms.Rank();
                 }
             }
          }, field_spec.typecode);
        }
        buffer_ptr += compound_type_size;
      }
    }
    
    hsize_t count[1] = {local_site_count};
    hid_t memspace_id = H5Screate_simple(1, count, NULL);
    hsize_t start[1] = {(hsize_t)comms.Scan(local_site_count, MPI_SUM) - local_site_count};
    H5_CHECK(H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5_CHECK(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE));
    H5_CHECK(H5Dwrite(dataset_id, compound_type_id, memspace_id, filespace_id, dxpl, buffer.data()));
    
    H5_CHECK(H5Pclose(dxpl));
    H5_CHECK(H5Sclose(memspace_id));
    H5_CHECK(H5Dclose(dataset_id));
    H5_CHECK(H5Sclose(filespace_id));
    H5_CHECK(H5Gclose(group_id));
    
    if (std::holds_alternative<single_timestep_files>(outputSpec.ts_mode) || timestepNumber >= totalSteps || (timestepNumber + outputSpec.frequency > totalSteps)) {
      if (file_id >= 0) {
        WriteXDMFFile();
        H5Fclose(file_id);
        file_id = -1; // Reset for next single file
      }
    }
  }
}
#endif
