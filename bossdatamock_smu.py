
def smufiles_noAP_binsplitted(catnamelist=catnamelist, catname2list=catname2list, RSDstrlist = ['noRSD', 'RSD'], 
                                  totbin=3, ibinlist = None, skipdata=False, onlydata=False):
    smufilelist = []
    if ibinlist == None: ibinlist = range(totbin)
    for catname in catnamelist:
        for ibin in ibinlist:
            if not skipdata:
                smufilelist.append(Tpcfrltfilename(binsplittedfilename(datafile(catname),  ibin+1, totbin)))
            if onlydata:
                continue
            for catname2 in catname2list:
                nummock = catinfo_nummock(catname, catname2)
                for imock in range(nummock):
                    for RSDstr in RSDstrlist:
                        smufilelist.append(Tpcfrltfilename(
                            binsplittedfilename(mockfile(catname, catname2, imock, RSDstr),  ibin+1, totbin)))
    return smufilelist

def smu__intxi__plot_all_binsplited_files(smu__intxi__settings=smu__intxi__settings_std, prefix='std-'):
    import bossdatamock; execfile(bossdatamock.pyfile)
    intxifilelist  = smu__intxi_calcwrite_list( 
                            smufilelist = smufiles_noAP_binsplitted(onlydata=True),  
                            smu__intxi__settings__list=[smu__intxi__settings],
                            just_get_intxifilelist=True
                            )
    time1 = time.time()
    intxis_list = [smu__intxi_quickload(nowfile) for nowfile in intxifilelist]
    time2 = time.time()
    mus = get_mid_array1d(np.linspace(0.01, 1.00, 41))
    fig, ax = figax(title = 'All data')
    for i in range(len(intxis_list)):
        ax.plot(mus, intxis_list[i])
    plt.show()
    fig.savefig(prefix+'Data.png', format='png')
    ifig = 1
    for catname in catnamelist:
        for nowcatname2list in [['HR3'], HR4LCcatname2list, ['HR4PSB']]:
            for RSDstr in RSDstrlist:
                for ibin in range(3):
                    intxifilelist  = smu__intxi_calcwrite_list( 
                            smufilelist = smufiles_noAP_binsplitted([catname], nowcatname2list, [RSDstr], 
                                                                    ibinlist=[ibin], skipdata=True),  
                            smu__intxi__settings__list=[smu__intxi__settings],
                            just_get_intxifilelist=True
                            )
                    time1 = time.time()
                    intxis_list = [smu__intxi_quickload(nowfile) for nowfile in intxifilelist]
                    time2 = time.time()
                    print time2 - time1
                    mus = get_mid_array1d(np.linspace(0.01, 1.00, 41))
                    figname = str_merge([catname, str_merge(nowcatname2list, ''), RSDstr, 'bin-'+str(ibin+1)], '-')
                    fig, ax = figax(title = figname)
                    for i in range(len(intxis_list)):
                        ax.plot(mus, intxis_list[i])
                    ax.set_ylim(5,28)
                    plt.show()
                    fig.savefig(prefix+figname+'.png', format='png')
                    ifig += 1
        print 'Done : ', ifig, 'figures.'


