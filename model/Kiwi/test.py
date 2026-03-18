import kiwi
print(kiwi.__file__)
kiwi.plot(
    gsn="gsn_metabolite_to_metabolite.txt",
    gss="GSS_Hacat_AB_metabo.txt",
    gls="GLS_Hacat_AB_metabo.txt",
    gsc="GSC_Hacat_AB_metabo.txt",
    pc=0.05, adj=True,
    hmt= 'values',
    gml='Hacat_AB_metabo.graphml',
    nwf='Hacat_AB_net_metabo.pdf',
    hmf='Hacat_AB_heatmap_metabo.pdf',
    lmp= True,
    saveHeatmapMatrix=True,
    heatmapMatrixFile='Hacat_AB_metabo_heatmap.csv'
)

kiwi.plot(
    gsn="gsn_pathway.txt",
    gss="GSS_Hacat_AB_pathway.txt",
    gls="GLS_Hacat_AB_pathway.txt",
    gsc="GSC_Hacat_AB_pathway.txt",
    pc=0.05, adj=True,
    hmt= 'values',
    gml='Hacat_AB_pathway.graphml',
    nwf='Hacat_AB_net_pathway.pdf',
    hmf='Hacat_AB_heatmap_pathway.pdf',
    lmp= True,
    saveHeatmapMatrix=True,
    heatmapMatrixFile='Hacat_AB_pathway_heatmap.csv'
)




