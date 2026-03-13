const { DOMParser } = require('@xmldom/xmldom');

const xml = `<?xml version="1.0" encoding="utf-8"?>
<Image TileSize="256" Overlap="1" Format="jpeg" xmlns="http://schemas.microsoft.com/deepzoom/2008">
  <Size Width="1000" Height="2000"/>
</Image>`;

try {
  const doc = new DOMParser().parseFromString(xml, "application/xml");
  const sizeElements = doc.getElementsByTagName("Size");
  
  if (sizeElements && sizeElements.length > 0) {
    const size = sizeElements[0];
    if (size.hasAttribute("Width")) {
      console.log("Width:", size.getAttribute("Width"));
    } else {
      console.warn("Size 標籤缺少 Width 屬性");
    }
  } else {
    console.warn("未找到 Size 標籤");
  }
} catch (error) {
  console.error("XML 解析失敗:", error);
}
