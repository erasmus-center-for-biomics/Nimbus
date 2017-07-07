
import collections
import naivetrack.utils
import naivetrack.models


FilterCriterion = collections.namedtuple(
    "FilterEntry", ["value", "quality_label", "depth_label"])


class SummarizeEntries(object):

    def __init__(self):
        """

        :return:
        """
        self.factors = [
            naivetrack.models.AlleleFactors["Sample"],
            naivetrack.models.AlleleFactors["Allele"]
            ]
        self.filters = []

    def add_filter(self, factor=None, equals="", quality_label=None, depth_label=None):
        """
        :param factor:
        :param quality_label:
        :param depth_label:

        :return:
        """
        assert factor is not None
        assert quality_label is not None
        assert depth_label is not None

        if len(self.factors) == 2:
            self.factors.append(int(factor))
        self.filters.append(
            FilterCriterion(str(equals), quality_label, depth_label)
        )

    def __call__(self, trackentry=None):
        """
        Determine the frequencies for a track entry

        :param trackentry:

        :return: a dict with the summaries added
        """
        assert trackentry is not None

        #
        sequences = set([a.sequence for a in trackentry.alleles])
        summaries = naivetrack.utils.aggregate_allele_information(
            trackentry, attributes=self.factors)

        rval = {}
        for sample in summaries.keys():

            # get the current summary
            summary = summaries[sample]

            # determine each data entry
            for sequence in sequences:

                # make sure to add a dict on the first sample
                if sequence not in rval.keys():
                    rval[sequence] = {}

                # add the entry
                if sample not in rval[sequence].keys():
                    rval[sequence][sample] = {}

                qval = 0
                dval = 0
                if len(self.factors) == 3 and sequence in summary.keys():
                    for x in summary[sequence].keys():
                        qval += summary[sequence][x]['__qual__']
                        dval += summary[sequence][x]['__depth__']
                elif sequence in summary.keys():
                    qval = summary[sequence]['__qual__']
                    dval = summary[sequence]['__depth__']

                rval[sequence][sample]['__quality__'] = qval
                rval[sequence][sample]['__depth__'] = dval

                if self.factors > 2 and sequence in summary.keys():
                    for i in xrange(0, len(self.filters)):
                        qval = 0
                        dval = 0
                        for x in summary[sequence].keys():
                            if x == self.filters[i].value:
                                qval += summary[sequence][x]['__qual__']
                                dval += summary[sequence][x]['__depth__']

                        rval[sequence][sample][self.filters[i].quality_label] = qval
                        rval[sequence][sample][self.filters[i].depth_label] = dval

        # return the summary data
        return rval
