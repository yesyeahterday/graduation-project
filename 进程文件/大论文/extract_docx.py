# -*- coding: utf-8 -*-
"""Extract text from docx for analysis"""
import zipfile
import xml.etree.ElementTree as ET
import re
import os

docx_path = r"c:\Users\23992\OneDrive\桌面\heyunlu\毕业设计\进程文件\大论文\复杂城市作战任务规划问题v1.docx"
output_path = r"c:\Users\23992\OneDrive\桌面\heyunlu\毕业设计\进程文件\大论文\extracted_ch4.txt"

# Extract text from docx (docx is a zip file)
with zipfile.ZipFile(docx_path, 'r') as z:
    xml_content = z.read('word/document.xml')
    
root = ET.fromstring(xml_content)
# Namespace for Word documents
ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}

paragraphs = []
for para in root.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}p'):
    texts = []
    for t in para.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t'):
        if t.text:
            texts.append(t.text)
        if t.tail:
            texts.append(t.tail)
    para_text = ''.join(texts).strip()
    if para_text:
        paragraphs.append(para_text)

# Also get table content
for tbl in root.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}tbl'):
    for row in tbl.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}tr'):
        row_texts = []
        for cell in row.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}tc'):
            cell_texts = []
            for t in cell.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t'):
                if t.text:
                    cell_texts.append(t.text)
                if t.tail:
                    cell_texts.append(t.tail)
            row_texts.append(''.join(cell_texts).strip())
        if any(row_texts):
            paragraphs.append(' | '.join(row_texts))

# Find Chapter 4 content
full_text = '\n'.join(paragraphs)
# Look for 第四章
ch4_start = full_text.find('第四章')
if ch4_start == -1:
    ch4_start = full_text.find('4.')
    
# Get content from Chapter 4 to end or next chapter (get more content for tables)
ch4_content = full_text[ch4_start:ch4_start+25000] if ch4_start >= 0 else full_text[:25000]

with open(output_path, 'w', encoding='utf-8') as f:
    f.write(ch4_content)

# Also save full doc for reference
full_path = output_path.replace('extracted_ch4.txt', 'extracted_full.txt')
with open(full_path, 'w', encoding='utf-8') as f:
    f.write(full_text[:30000])  # First 30000 chars

print("Extracted to:", output_path)
print("Length:", len(ch4_content))
