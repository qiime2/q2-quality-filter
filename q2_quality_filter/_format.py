import qiime2.plugin.model as model


class QualityFilterStatsFmt(model.TextFileFormat):
    def sniff(self):
        line = open(str(self)).readline()
        hdr = line.strip().split(',')
        expected = ['sample-id', 'total-input-reads', 
                    'reads-truncated', 
                    'reads-too-short-after-truncation', 
                    'reads-exceeding-maximum-ambiguous-bases']
        return hdr == expected


QualityFilterStatsDirFmt = model.SingleFileDirectoryFormat(
    'QualityFilterStatsDirFmt', 'stats.csv', QualityFilterStatsFmt)
