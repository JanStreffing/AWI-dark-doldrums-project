#!/bin/bash
# Self-contained, session-independent chain: finish 2025, then all of 2024, then all of 2023.
cd /work/ab0246/a270092/software/dark-doldrums-project/solar

echo "=== $(date) starting 2025 remaining ===" >> logs/download_chain.log
./scripts/download/download_era5.sh 202501 202502 202503 202504 202505 202506 202507 202509 202510 >> logs/download_chain.log 2>&1

echo "=== $(date) starting 2024 (full year) ===" >> logs/download_chain.log
./scripts/download/download_era5.sh 202401 202402 202403 202404 202405 202406 202407 202408 202409 202410 202411 202412 >> logs/download_chain.log 2>&1

echo "=== $(date) starting 2023 (full year) ===" >> logs/download_chain.log
./scripts/download/download_era5.sh 202301 202302 202303 202304 202305 202306 202307 202308 202309 202310 202311 202312 >> logs/download_chain.log 2>&1

echo "=== $(date) CHAIN COMPLETE ===" >> logs/download_chain.log
