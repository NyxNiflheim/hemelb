#!/bin/bash


# --- 配置区 ---
HEMELB_EXECUTABLE="./hemelb/build/hemelb"
H5_CONFIG_FILE="./hemelbTest-data/tests/gmy/bifurcation-40um/configH5.xml"
XTR_CONFIG_FILE="./hemelbTest-data/tests/gmy/bifurcation-40um/configXTR.xml"
# 脚本日志文件的输出目录
BENCHMARK_LOG_DIR="./benchmark_run" 
# HemeLB生成results和report.txt的目录
SIMULATION_CASE_DIR="./hemelbTest-data/tests/gmy/bifurcation-40um"
RESULTS_DIR="${SIMULATION_CASE_DIR}/results"

# 要测试的核心数列表
declare -a CORE_COUNTS=(4)


echo "--- Starting HemeLB I/O Benchmark ---"

mkdir -p ${BENCHMARK_LOG_DIR}

for format in "XTR" "Hdf5"; do
  for cores in "${CORE_COUNTS[@]}"; do

    # 根据格式选择正确的配置文件
    if [ "$format" == "XTR" ]; then
      CONFIG_TO_USE=${XTR_CONFIG_FILE}
      OUTPUT_FILENAME="${RESULTS_DIR}/Extracted/flow_snapshot.dat"
    else
      CONFIG_TO_USE=${H5_CONFIG_FILE}
      OUTPUT_FILENAME="${RESULTS_DIR}/Extracted/flow_snapshot.h5"
    fi
    
    # # 从配置文件中动态读取硬编码的输出文件名
    # OUTPUT_FILENAME="${RESULTS_DIR}/Extracted/flow_snapshot.h5"

    # 运行前清理旧的results目录
    rm -rf ${RESULTS_DIR}
    
    # 运行模拟
    mpirun -np ${cores} ${HEMELB_EXECUTABLE} -in ${CONFIG_TO_USE} > "${BENCHMARK_LOG_DIR}/run_${format}_${cores}.log" 2>&1
    
    # report
    REPORT_FILE="${RESULTS_DIR}/report.txt"
    if [ -f "$REPORT_FILE" ]; then
        line=$(grep "Extraction writing" "$REPORT_FILE")
        time1=$(echo "$line" | awk '{print $3}')
        time2=$(echo "$line" | awk '{print $4}')
        time3=$(echo "$line" | awk '{print $5}')
        time4=$(echo "$line" | awk '{print $6}')
        echo "${format},process number ${cores}"
        echo "  Extraction writing ${time1},${time2},${time3},${time4}"
    else
        echo "${format},${cores},ERROR_NO_REPORT"
    fi

    if [ -f "$OUTPUT_FILENAME" ]; then
        FILE_SIZE_BYTES=$(ls -l ${OUTPUT_FILENAME} | awk '{print $5}')
        FILE_SIZE_MB=$(echo "scale=2; ${FILE_SIZE_BYTES} / 1024 / 1024" | bc)
        echo "  File Size ${FILE_SIZE_MB} MB"
    else
        echo "  ERROR_NO_FILE"
    fi

  done
done

echo "--- Benchmark Finished ---"

