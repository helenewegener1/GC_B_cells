for f in ../out/*/idmap.txt; do
  dir=$(basename "$(dirname "$f")")
  cp "$f" "./idmap_${dir}.txt"
done
