$MGT_HOME/test_data/two_vir_fam

        use_sign=False
        normalization=FULL_NORMALIZATION
        feats_train,feats_test = seqToWordFeatures(dataTrain,dataTest,order,gap)
        kernel=kernelClass(
            feats_train, feats_train, use_sign, normalization)
