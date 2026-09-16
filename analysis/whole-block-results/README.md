# 250 公尺選取完整街廓版

新規則：市民大道原 250 m 雙側範圍碰到的完整街廓全部保留，再與原範圍取聯集以保留道路。面積 7.396 km²；新增 1.419 km²。

街廓以現有 OpenStreetMap 地面道路網形成的封閉面推定，共 932 面；不使用建築量體，也不讓高架橋分割街廓。這是道路中心線街廓近似值，並非地籍或實測路緣。缺失道路不以猜測補線；本資料無法保證尚未形成封閉面的街廓已完整辨認。

- [完整向量 SVG](whole-block-vegetation.svg.gz)：綠色植被、土黃待核實、紫灰新增未判讀範圍、淺灰原分類其他。SVG 只有向量路徑，沒有嵌入點陣影像。
- [街廓與原 250 m 邊界 SVG](block-boundary.svg)
- [新研究邊界 GeoJSON](study-boundary.geojson)
- [完整街廓 GeoJSON](retained-blocks.geojson)
- [分類向量 GeoJSON（gzip）](vegetation.geojson.gz)：本機保存；GitHub 提供 SVG 壓縮副本（解壓後使用）。

EPSG:3826 以公尺等比例製作 SVG；GeoJSON 為 WGS84。新增範圍尚未做植被分類，因此舊 10.36% 不適用此版新邊界。本次完成邊界與向量調整，沒有假設新增面積為零綠覆。
