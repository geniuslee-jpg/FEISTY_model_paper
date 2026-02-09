#!/bin/bash
# =============================================================
# 00_preprocess_cdo.sh
# NEMO-FEISTY 출력에서 필요한 변수만 추출 + 10년 시간 자르기
#
# 사용법:
#   cd /data01/labdisk/sungjin/NEMO-FEISTY
#   bash 00_preprocess_cdo.sh
#
# 출력: ./subset/ 디렉토리에 전처리된 파일 생성
# =============================================================

set -e

# --- 설정 ---
INDIR="."
OUTDIR="./subset"
mkdir -p ${OUTDIR}

echo "=============================================="
echo " NEMO-FEISTY CDO 전처리"
echo " 입력: ${INDIR}"
echo " 출력: ${OUTDIR}"
echo "=============================================="

# =============================================================
# 1. grid_T_static (시간 없음, 변수만 추출)
#    필요: depth (해저 깊이)
# =============================================================
echo ""
echo "[1/4] grid_T_static → depth만 추출..."
IN="${INDIR}/ORCA2_1y_00010101_00501231_grid_T_static.nc"
OUT="${OUTDIR}/static_depth.nc"

cdo selvar,depth ${IN} ${OUT}
echo "  → $(du -h ${OUT} | cut -f1)  ${OUT}"

# =============================================================
# 2. grid_T 5d → thetao, sbt / 10년 (timestep 1-730)
#    noleap calendar, 5일 간격: 73 steps/yr × 10yr = 730 steps
#    필요: thetao (3D 온도), sbt (해저 온도)
# =============================================================
echo ""
echo "[2/4] grid_T 5d → thetao + sbt, 10년 추출..."
IN="${INDIR}/ORCA2_5d_00010101_00501231_grid_T.nc"

# Step 2a: 변수 추출
TMP1="${OUTDIR}/_tmp_gridT_vars.nc"
cdo selvar,thetao,sbt ${IN} ${TMP1}
echo "  변수 추출 완료 (thetao, sbt)"

# Step 2b: 10년 시간 추출
OUT="${OUTDIR}/grid_T_10yr.nc"
cdo seltimestep,1/730 ${TMP1} ${OUT}
rm -f ${TMP1}
echo "  → $(du -h ${OUT} | cut -f1)  ${OUT}"

# =============================================================
# 3. diad_T 1y → Heup, EPC100, GRAZ1, GRAZ2 / 10년 (timestep 1-10)
#    필요: Heup (유광층), EPC100 (export), GRAZ1/GRAZ2 (grazing)
# =============================================================
echo ""
echo "[3/4] diad_T 1y → Heup + EPC100 + GRAZ1 + GRAZ2, 10년..."
IN="${INDIR}/ORCA2_1y_00010101_00501231_diad_T.nc"

# 변수 추출 + 시간 추출 (파이프)
OUT="${OUTDIR}/diad_T_1y_10yr.nc"
cdo seltimestep,1/10 -selvar,Heup,EPC100,GRAZ1,GRAZ2 ${IN} ${OUT}
echo "  → $(du -h ${OUT} | cut -f1)  ${OUT}"

# =============================================================
# 4. diad_T 1m → 같은 변수 / 10년 (timestep 1-120, 12개월×10년)
#    (월별 데이터도 백업으로 보관)
# =============================================================
echo ""
echo "[4/4] diad_T 1m → Heup + EPC100 + GRAZ1 + GRAZ2, 10년..."
IN="${INDIR}/ORCA2_1m_00010101_00501231_diad_T.nc"

OUT="${OUTDIR}/diad_T_1m_10yr.nc"
cdo seltimestep,1/120 -selvar,Heup,EPC100,GRAZ1,GRAZ2 ${IN} ${OUT}
echo "  → $(du -h ${OUT} | cut -f1)  ${OUT}"

# =============================================================
# 요약
# =============================================================
echo ""
echo "=============================================="
echo " 전처리 완료!"
echo "=============================================="
echo ""
echo " 생성된 파일:"
ls -lh ${OUTDIR}/
echo ""
echo " 원본 대비 크기:"
echo "   원본 합계: $(du -sh ${INDIR}/*.nc | tail -1 | cut -f1)"
echo "   subset 합계: $(du -sh ${OUTDIR}/ | cut -f1)"
echo ""
echo " 다음 단계:"
echo "   python 02_extract_feisty_input.py ${OUTDIR}"
