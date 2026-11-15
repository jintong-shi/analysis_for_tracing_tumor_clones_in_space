#!/usr/bin/env bash
set -euo pipefail

# ==== 统一根目录（你要的地方）====
BASE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$BASE"/{cache,tmp,conda,pkgs,pip,mpl,rcache,mamba}

# 常见缓存/临时目录全部指过去
export XDG_CACHE_HOME="$BASE/cache"
export TMPDIR="$BASE/tmp"
export PIP_CACHE_DIR="$BASE/pip"
export MPLCONFIGDIR="$BASE/mpl"
export R_USER_CACHE_DIR="$BASE/rcache"
export CONDA_PKGS_DIRS="$BASE/pkgs"

# micromamba 的根（环境和包缓存都会放这里）
export MAMBA_ROOT_PREFIX="$BASE/mamba"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

# ==== 安装 micromamba（若无）====
if ! command -v micromamba >/dev/null 2>&1; then
  echo "[bootstrap] Installing micromamba (linux-x86_64) into $MAMBA_ROOT_PREFIX ..."
  mkdir -p "$MAMBA_ROOT_PREFIX"
  curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "$MAMBA_ROOT_PREFIX" bin/micromamba
fi

# 启用 shell hook（脚本内需，脚本外你仍需手动 activate）
eval "$(micromamba shell hook -s bash)"

# ==== 创建运行 snakemake 的小环境 ====
ENV_NAME="snakemake-run"
if ! micromamba env list | grep -qE "^\s*$ENV_NAME\s"; then
  micromamba create -y -n "$ENV_NAME" -c conda-forge -c bioconda snakemake git
fi

# （可选）脚本内临时激活以创建目录等
micromamba activate "$ENV_NAME"

# 预创建项目目录
mkdir -p config envs R patches software envflags

# 补丁脚本赋权
[ -f patches/copyvae_postinstall.sh ] && chmod +x patches/copyvae_postinstall.sh || true

echo "[bootstrap] Done.

Next steps:
  1) micromamba activate $ENV_NAME
  2) snakemake -j 8 --use-conda --conda-frontend micromamba \\
       --conda-prefix \"$BASE/conda\"
"