#!/usr/bin/env bash
set -euo pipefail

cd software/copyVAE

# 已打过补丁就跳过
if grep -q "GENE_MAP_CSV" copyvae/pipeline.py; then
  echo "[patch] copyVAE already patched for offline gene_map. Skip."
  exit 0
fi

PATCH_REL="../../patches/copyvae_pipeline_offline_genemap.patch"

echo "[patch] Applying copyVAE offline gene_map patch..."
# 优先使用 git apply；失败则回退到 patch 命令
if git apply --check "${PATCH_REL}" >/dev/null 2>&1; then
  git apply --whitespace=fix "${PATCH_REL}"
elif command -v patch >/dev/null 2>&1; then
  patch -p1 < "${PATCH_REL}"
else
  echo "[patch] ERROR: neither 'git apply' nor 'patch' worked." >&2
  exit 1
fi

# 开发模式安装已经是 -e，不必重新 pip install；如需可解除注释：
# pip install -e . --no-deps

echo "[patch] Patch applied."