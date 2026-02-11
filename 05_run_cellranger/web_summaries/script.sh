# More efficient: directly target the expected structure
for sample_dir in ../out_v9/res_* ../out_v8/res_*; do
    if [ -d "$sample_dir" ]; then
        # Look for web_summary.html in the expected location
        html_file="$sample_dir/outs/per_sample_outs/$(basename $sample_dir)/web_summary.html"
        
        if [ -f "$html_file" ]; then
            sample=$(basename "$sample_dir")
            cp "$html_file" "web_summaries/${sample}_web_summary.html"
            echo "Copied: $sample"
        fi
    fi
done
