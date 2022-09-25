## 工具源代码网址

| Methods       | Type                     | Language | source                                               | achieve      |
| ------------- | ------------------------ | -------- | ---------------------------------------------------- | ------------ |
| SPCS          | Diffusion                | R        | https://github.com/Usos/SPCS                         | inline       |
| SpotClean     | Decontamination          | R        | https://github.com/zijianni/SpotClean                | Bioconductor |
| stlearn       | Clustering               | Python   | https://github.com/BiomedicalMachineLearning/stLearn | pip          |
| SpaGCN        | Clustering, Finding-HVGS | Python   | https://github.com/jianhuupenn/SpaGCN                | pip          |
| BayesSpace    | Clustering               | R        | https://github.com/edward130603/BayesSpace           | Bioconductor |
| SomDE         | Clustering, Finding-HVGS | Python   | https://github.com/WhirlFirst/somde                  | pip          |
| SpatialDE     | Finding-HVGs             | Python   | https://github.com/PMBio/SpatialDE                   | pip          |
| SPARK         | Finding-HVGs             | R        | https://github.com/xzhoulab/SPARK                    | github       |
| trendsceek    | Finding-HVGs             | R        | https://github.com/edsgard/trendsceek                | github       |
| SPOTlight     | Deconvolution            | R        | https://github.com/MarcElosua/SPOTlight              | github       |
| spacexr       | Deconvolution            | R        | https://github.com/dmcable/spacexr                   | github       |
| stereoScope   | Deconvolution            | Python   | https://github.com/almaan/stereoscope                | inline       |
| Cell2Location | Deconvolution            | Python   | https://github.com/BayraktarLab/cell2location        | pip          |
| Giotto        | Pipline                  | R        | https://github.com/RubD/Giotto/                      | github       |
| squidpy       | Pipline                  | Python   | https://github.com/theislab/squidpy                  | pip          |
| NNLM          | tool                     | R        | https://github.com/linxihui/NNLM                     | github       |
| SpaTalk       | Cell-Cell interaction    | R        | https://github.com/ZJUFanLab/SpaTalk                 | github       |
| iTalk         | Cell-Cell interaction    | R        | https://github.com/Coolgenome/iTALK                  | github       |
| CellChat      | Cell-Cell interaction    | R        | https://github.com/sqjin/CellChat                    | Bioconductor |
| CellCall      | Cell-Cell interaction    | R        | https://github.com/ShellyCoder/cellcall              | github       |

## 信号通路全解

### FGF

成纤维细胞生长因子（fibroblast growth factor, FGF）有几种异构体，在[动脉硬化](https://baike.baidu.com/item/动脉硬化/1521831?fromModule=lemma_inlink)灶中起作用的主要是[bFGF](https://baike.baidu.com/item/bFGF/5573258?fromModule=lemma_inlink)（basic fibroblast growth factor），bFGF可以由[内皮细胞](https://baike.baidu.com/item/内皮细胞/5911806?fromModule=lemma_inlink)、[平滑肌细胞](https://baike.baidu.com/item/平滑肌细胞/1333038?fromModule=lemma_inlink)、[巨噬细胞](https://baike.baidu.com/item/巨噬细胞/245209?fromModule=lemma_inlink)分泌。它的作用是促进内皮细胞的游走和平滑肌细胞的增殖，不能使平滑肌细胞游走。能够促进新血管形成，修复损害的内皮细胞。FGF被认为是病灶形成促进因子，但从修复角度看它也有有利的一面。

####  FGF信号通路

FGFs家族成员包括FGF1-23,根据FGFs的作用机制可将其分为3类：FGF11-14为非分泌信号且不依赖于FGF受体(fibroblast growth factor receptors,  FGFRs)而发挥作用,可以影响神经元的电兴奋性,被称为胞内FGFs;其余两类FGFs皆依赖于FGFRs,通过局部信号的扩散而作用于附近的靶细胞,主要作为胚胎发育中组织形态和器官形成的旁分泌因子(包括FGF1/2/5, FGF3/4/6, FGF7/10/22,  FGF8/17/18和FGF9/16/20亚家族)和多种代谢过程中的内分泌信号的内分泌因子(FGF15/19/21/23亚家族)[[6](javascript:;)]。在FGFs信号通路中,FGFs与细胞膜表面FGFRs结合标志着FGF信号通路的激活。目前有4种包含细胞内酪氨酸激酶结构域的FGFRs  (FGFR1~4)。许多研究发现FGFRs细胞内区域的酪氨酸磷酸化可激活细胞内关键的信号通路,包括促分裂素原活化蛋白激酶(mitogenactivated protein kinase, MAPK)、磷脂酰肌醇3-激酶-Akt(phosphoinositide 3-kinase-AKT,  PI3K-Akt)、磷脂酶Cγ/蛋白激酶Cα(phospholipase Cγ/proteinkinase Cα,  PLCγ/PKCα)和信号转导子与转录激活子(signal transducer and activator of transcription, STAT)信号通路[[7](javascript:;)]。FGF信号通路在肝脏、肾脏、脑部和骨骼等组织器官的发育以及生理和病理过程中均发挥重要作用,内耳的发育也与FGF信号通路密切相关。

### TGF-β

体内组织中的细胞增殖，胚胎发生、分化和细胞死亡过程中细胞的特定命运都受到多种细胞与细胞之间信号的控制，这种控制一旦发生异常，将会带来非常严重的后果。这些调控信号中最突出的是TGF-Beta超家族，该家族包含大量不同的多肽形态发生因子，包括TGF-Beta本身以及BMP（骨形态发生蛋白）和GDF（增长和分化因子）（参考文献1）。TGF-Beta家族的成员会在不同的时间点并以组织特异性的形式表达，因此在机体中大多数组织的发育、稳态平衡和修复中起重要作用。所有免疫细胞，包括B细胞、T细胞和树突状细胞，以及巨噬细胞，都分泌TGF-Beta，而TGF-Beta又通过其他细胞因子负调控免疫细胞的增殖、分化和激活。 **因此，TGF-Beta是一种有效的免疫抑制剂，TGF-β信号的紊乱还与自身免疫、炎症和癌症有关**（参考文献2）。

\1. Ripamonti U, Crooks J, Matsaba T,  TaskerJ. Induction of endochondral bone formation by recombinant  humantransforming growth factor-beta2 in the baboon (Papio ursinus).  GrowthFactors. 2000; 17(4): 269-85. PubMed ID: 10801076

\2. Moustakas A, Pardali K, Gaal A,  HeldinCH. Mechanisms of TGF-beta signaling in regulation of cell growth  anddifferentiation. Immunol Lett. 2002 Jun 3; 82(1-2): 85-91. PubMed ID: 12008039

### CXCL

研究发现CXCL1/CXCR2通路介导肿瘤间质微环境，在胰腺癌模型中阻断该通路可以延长生存期；在胃癌中该通路与肿瘤转移、分期和预后密切相关，调控淋巴管内皮细胞增生和功能。我们推测胃癌细胞通过CXCL1/CXCR2通路招募和改造间质细胞，促进自身生长和播散。本课题将利用临床标本预筛CXCL1/CXCR2通路调控的间质细胞；应用IHC、qPCR和ELISA体外研究CXCL1表达对细胞表型和功能的作用，并研究作用机制；通过皮下移植瘤模型研究间质细胞CXCR2表达对胃癌细胞成瘤性和瘤特征的影响；揭示胃癌细胞分泌CXCL1,依赖与间质细胞表面CXCR2结合来调控间质细胞表型和功能，促进胃癌增生和转移；最后，利用胃癌移植瘤模型，研究该通路阻断剂联合细胞毒药物抑制胃癌增生和转移的效果。本研究旨在开辟胃癌治疗新靶点，改善胃癌患者预后。

CXCL12/CXCR4



### PI3K/AKT/mTOR

https://zhuanlan.zhihu.com/p/115590521

