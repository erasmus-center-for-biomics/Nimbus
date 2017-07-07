import table.Location
import table

class LocationReader(table.HeaderTable):
    """
    LocationReader adds location data to a header.

    Wraps a table object associate a
    Location to the values.
    """

    def __init__(self, *args, **kwargs):
        """__init__ creates a new LocationParser."""
        super(LocationReader, self).__init__(*args, **kwargs)
        self.location_indexes = []

        # add the chromosome list
        self.chromosome_list = []
        if "chromosome_list" in kwargs.keys():
            self.chromosome_list = kwargs["chromosome_list"]

        # labels{chromosome, start, end}
        self.labels = [None, None, None]
        if "label_chr" in kwargs.keys():
            self.labels[0] = kwargs["label_chr"]
        if "label_start" in kwargs.keys():
            self.labels[1] = kwargs["label_start"]
        if "label_end" in kwargs.keys():
            self.labels[2] = kwargs["label_end"]

        # make sure the labels have been added
        if self.labels[0] is None:
            raise ValueError("chromosome label not provided (label_chr)")
        if self.labels[1] is None:
            raise ValueError("chromosome label not provided (label_start)")
        if self.labels[2] is None:
            raise ValueError("chromosome label not provided (label_end)")

        # is the table one-based or zero-based (return data
        # is always zero based)
        self.one_based = True if "one_based" in kwargs.keys() else False

    def next(self):
        """Return the next element."""
        values = super(LocationReader, self).next()

        # add the labels if not previously defined
        if len(self.location_indexes) == 0:
            # add the labels
            for idx in xrange(len(self.labels)):
                if self.labels[idx] in self.header:
                    self.location_indexes.append(idx)
                else:
                    raise ValueError("label %s not found in header" % self.labels[idx])

        # parse the chromosome information
        try:
            chridx = self.chromosome_list.index(values[self.location_indexes[0]][1])
        except ValueError:
            self.chromosome_list.append(values[self.location_indexes[0]][1])
            chridx = len(self.chromosome_list) - 1

        # parse the location
        start = int(values[self.location_indexes[1]][1])
        end = int(values[self.location_indexes[2]][1])
        if self.one_based:
            start -= 1
            end -= 1

        # set the location
        location = table.Location.Loc(chridx, start, end)

        # return the location and the values
        return location, values


class LocationMultiAnnoReader(table.MultiAnno):
    """
    LocationMultiAnnoReader adds location data to a header.

    Wraps a table object associate a
    Location to the values.
    """

    def __init__(self, *args, **kwargs):
        """__init__ creates a new LocationMultiAnnoReader."""
        super(LocationMultiAnnoReader, self).__init__(*args, **kwargs)
        self.location_indexes = []

        # add the chromosome list
        self.chromosome_list = []
        if "chromosome_list" in kwargs.keys():
            self.chromosome_list = kwargs["chromosome_list"]

        # labels{chromosome, start, end}
        self.labels = [None, None, None]
        if "label_chr" in kwargs.keys():
            self.labels[0] = kwargs["label_chr"]
        if "label_start" in kwargs.keys():
            self.labels[1] = kwargs["label_start"]
        if "label_end" in kwargs.keys():
            self.labels[2] = kwargs["label_end"]

        # make sure the labels have been added
        if self.labels[0] is None:
            raise ValueError("chromosome label not provided (label_chr)")
        if self.labels[1] is None:
            raise ValueError("chromosome label not provided (label_start)")
        if self.labels[2] is None:
            raise ValueError("chromosome label not provided (label_end)")

        # is the table one-based or zero-based (return data
        # is always zero based)
        self.one_based = True if "one_based" in kwargs.keys() else False

    def next(self):
        """Return the next element."""
        values = super(LocationMultiAnnoReader, self).next()

        # add the labels if not previously defined
        if len(self.location_indexes) == 0:
            # add the labels
            for idx in xrange(len(self.labels)):
                if self.labels[idx] in self.header:
                    self.location_indexes.append(idx)
                else:
                    raise ValueError("label %s not found in header" % self.labels[idx])

        # parse the chromosome information
        try:
            chridx = self.chromosome_list.index(values[self.location_indexes[0]][1])
        except ValueError:
            self.chromosome_list.append(values[self.location_indexes[0]][1])
            chridx = len(self.chromosome_list) - 1

        # parse the location
        start = int(values[self.location_indexes[1]][1])
        end = int(values[self.location_indexes[2]][1])
        if self.one_based:
            start -= 1
            end -= 1

        # set the location
        location = table.Location.Loc(chridx, start, end)

        # return the location and the values
        return location, values

