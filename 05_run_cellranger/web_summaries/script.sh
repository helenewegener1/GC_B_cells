# More efficient: directly target the expected structure
for sample_dir in ../out_v9/res_*; do
   
    # Look for web_summary.html in the expected location
    html_file="$sample_dir/outs/per_sample_outs/$(basename $sample_dir)/web_summary.html"
    sample=$(basename "$sample_dir")
    cp "$html_file" "./${sample}_web_summary.html"
    
done
