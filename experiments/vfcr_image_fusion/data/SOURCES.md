# 数据来源

- BSDS300：论文明确指定的 Berkeley Segmentation Dataset 300。
- `385028.jpg`：论文图 5 中的房屋图像。
- `113016.jpg`：论文图 6 中的马图像。
- `boat.512.tiff`：USC-SIPI 的 512×512 灰度船图，与论文图 8 的原图一致。
- `t1_icbm_normal_1mm_pn0_rf0.rawb`：BrainWeb 的 T1、1 mm、0% 噪声、0% 灰度不均匀体数据。

下载地址：

- https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/
- https://sipi.usc.edu/database/
- https://brainweb.bic.mni.mcgill.ca/brainweb/

`download_data.ps1` 会下载原始压缩包并保留在 `data/raw` 中。
