import cmath
import matplotlib.pyplot as plt
import csv

CONST_VALID_REQUESTS = ['ACCELERATION', 'DISPLACEMENTS', 'MPCF', 'SPCF', 'ELEMENT FORCES', 'ELEMENT STRAINS']


def dispatch_parse(output, data_chunks):
    if output == 'MAGNITUDE-PHASE' or output == 'REAL-IMAGINARY':
        num = int(len(data_chunks) / 2)
        if len(data_chunks) % 2 != 0:
            raise ValueError('Wrong number of chunks!', 'Output: %s, num of chunks: %d' % (output, len(data_chunks)))
    else:
        num = len(data_chunks)

    if output == 'MAGNITUDE-PHASE':
        return [data_chunks[i]*cmath.exp(1j*data_chunks[i+num]*cmath.pi/180.0) for i in range(num)]
    elif output == 'REAL-IMAGINARY':
        return [data_chunks[i] + 1j*data_chunks[i+num] for i in range(num)]
    else:
        return [data_chunks[i] for i in range(num)]


class PchParser:
    def reset_current_frame(self):
        self.cur_data_chunks = []
        self.is_frequency_response = False
        self.output_sort = 0
        self.cur_subcase = 0
        self.cur_output = 0
        self.current_frequency = 0
        self.cur_entity_id = 0
        self.cur_entity_type_id = 0

    def __init__(self, filename):
        # define the dictionary
        self.parsed_data = {'FREQUENCY': {}, 'SUBCASES': set()}
        for request in CONST_VALID_REQUESTS:
            self.parsed_data[request] = {}

        # initiate current frame
        self.reset_current_frame()

        is_header = True

        # start reading
        with open(filename, 'r') as pch:
            # read only first 72 characters from the punch file
            for line in pch:
                line = line[0:72]

                # reset all variables
                if line.startswith('$TITLE   ='):
                    is_header = False
                    # insert the last frame remaining in memory
                    self.insert_current_frame()
                    # reset the frame
                    self.reset_current_frame()

                # skip everything before TITLE
                if is_header:
                    continue

                # parse the subcase
                if line.startswith('$SUBCASE ID =') or \
                        line.startswith('$RANDOM ID ='):
                    self.cur_subcase = int(line[13:].strip())
                    self.parsed_data['SUBCASES'].add(self.cur_subcase)

                # identify NASTRAN request
                if line.startswith('$DISPLACEMENTS'):
                    self.cur_request = 'DISPLACEMENTS'
                elif line.startswith('$ACCELERATION'):
                    self.cur_request = 'ACCELERATION'
                elif line.startswith('$MPCF'):
                    self.cur_request = 'MPCF'
                elif line.startswith('$SPCF'):
                    self.cur_request = 'SPCF'
                elif line.startswith('$ELEMENT FORCES'):
                    self.cur_request = 'ELEMENT FORCES'
                elif line.startswith('$ELEMENT STRAINS'):
                    self.cur_request = 'ELEMENT STRAINS'

                # identify output type
                if line.startswith('$REAL-IMAGINARY OUTPUT'):
                    self.cur_output = 'REAL-IMAGINARY'
                elif line.startswith('$MAGNITUDE-PHASE OUTPUT'):
                    self.cur_output = 'MAGNITUDE-PHASE'
                elif line.startswith('REAL OUTPUT'):
                    self.cur_output = 'REAL'

                # parse of frequency response results
                if line.find('IDENTIFIED BY FREQUENCY') != -1:
                    self.is_frequency_response = True
                    self.output_sort = 2
                elif line.find('$FREQUENCY =') != -1:
                    self.is_frequency_response = True
                    self.output_sort = 1

                # parse entity id
                if line.startswith('$POINT ID ='):
                    self.cur_entity_id = int(line[11:23].strip())
                elif line.startswith('$ELEMENT ID ='):
                    self.cur_entity_id = int(line[13:23].strip())
                elif line.startswith('$FREQUENCY = '):
                    self.current_frequency = float(line[12:28].strip())

                # parse element type
                if line.startswith('$ELEMENT TYPE ='):
                    self.cur_entity_type_id = int(line[15:27].strip())

                # ignore other comments
                if line.startswith('$'):
                    continue

                # check if everything ok
                self.validate()

                # start data parsing
                line = line.replace('G', ' ')
                if line.startswith('-CONT-'):
                    line = line.replace('-CONT-', '')
                    self.cur_data_chunks += [float(_) for _ in line.split()]
                else:
                    # insert the last frame
                    self.insert_current_frame()

                    # update the last frame with a new data
                    self.cur_data_chunks = [float(_) for _ in line.split()]

            # last block remaining in memory
            self.insert_current_frame()

    def validate(self):
        if self.cur_request not in CONST_VALID_REQUESTS:
            raise NotImplementedError("Request %s is not implemented", self.cur_request)

        if self.cur_request == 'ELEMENT FORCES' and self.cur_entity_type_id not in [12, 102]:
            raise NotImplementedError("Element forces parser is implemented only for CELAS2 and CBUSH elements!")

    def insert_current_frame(self):
        # last block remaining in memory
        if len(self.cur_data_chunks) > 0:
            # ensure that subcase is allocated in the dataset
            if self.cur_subcase not in self.parsed_data[self.cur_request]:
                self.parsed_data[self.cur_request][self.cur_subcase] = {}
                self.parsed_data['FREQUENCY'][self.cur_subcase] = {}

            values = dispatch_parse(self.cur_output, self.cur_data_chunks[1:])
            if self.is_frequency_response:
                # incremented by frequency, entity is given
                if self.output_sort == 2:
                    self.current_frequency = self.cur_data_chunks[0]
                # incremented by entity, frequency is given
                elif self.output_sort == 1:
                    self.cur_entity_id = int(self.cur_data_chunks[0])

                # insert frequency in the database
                if self.current_frequency not in self.parsed_data['FREQUENCY'][self.cur_subcase]:
                    self.parsed_data['FREQUENCY'][self.cur_subcase][self.current_frequency] = \
                        len(self.parsed_data['FREQUENCY'][self.cur_subcase])

                # ensure that dictionary for the entity exists
                if self.cur_entity_id not in self.parsed_data[self.cur_request][self.cur_subcase]:
                    self.parsed_data[self.cur_request][self.cur_subcase][self.cur_entity_id] = []

                self.parsed_data[self.cur_request][self.cur_subcase][self.cur_entity_id].append(values)
            else:
                self.cur_entity_id = int(self.cur_data_chunks[0])
                self.parsed_data[self.cur_request][self.cur_subcase][self.cur_entity_id] = values

    def health_check(self):
        frequency_steps = []
        for subcase in self.parsed_data['SUBCASES']:
            frequency_steps.append(len(self.parsed_data['FREQUENCY'][subcase]))
        assert min(frequency_steps) == max(frequency_steps)

    def get_subcases(self):
        return sorted(self.parsed_data['SUBCASES'])

    def __get_data_per_request(self, request, subcase):
        self.health_check()
        if subcase in self.parsed_data[request]:
            return self.parsed_data[request][subcase]
        else:
            raise KeyError('%s data for subase %s is not found' % (request, subcase))

    def get_accelerations(self, subcase, entityID=None, component=None):
        if (entityID != component) and (entityID == None or component == None):
            raise KeyError("Need both Entity ID and direction component")

        elif entityID == None:
            return self.__get_data_per_request('ACCELERATION', subcase)

        elif component == 'tx':
            tx_data = []
            for each in self.get_accelerations(subcase)[entityID]:
                tx_data.append(each[0])
            return tx_data

        elif component == 'ty':
            ty_data = []
            for each in self.get_accelerations(subcase)[entityID]:
                ty_data.append(each[1])
            return ty_data

        elif component == 'tz':
            tz_data = []
            for each in self.get_accelerations(subcase)[entityID]:
                tz_data.append(each[2])
            return tz_data

        elif component == 'rx':
            rx_data = []
            for each in self.get_accelerations(subcase)[entityID]:
                rx_data.append(each[3])
            return rx_data

        elif component == 'ry':
            ry_data = []
            for each in self.get_accelerations(subcase)[entityID]:
                ry_data.append(each[4])
            return ry_data

        elif component == 'rz':
            rz_data = []
            for each in self.get_accelerations(subcase)[entityID]:
                rz_data.append(each[5])
            return rz_data

    def get_displacements(self, subcase):
        return self.__get_data_per_request('DISPLACEMENTS', subcase)

    def get_mpcf(self, subcase):
        return self.__get_data_per_request('MPCF', subcase)

    def get_spcf(self, subcase):
        return self.__get_data_per_request('SPCF', subcase)

    def get_forces(self, subcase):
        return self.__get_data_per_request('ELEMENT FORCES', subcase)

    def get_frequencies(self, subcase):
        return sorted(self.parsed_data['FREQUENCY'][subcase])


class SimplePch:
    """Parser for a simple punch file with no headers. Useful for random
    vibration punch output.
    """

    def __init__(self, filename):
        # initiate dictionary
        self.entitylist = []
        self.data = {}

        # start reading
        with open(filename, 'r') as pch:
            # read only first 72 characters from the punch file
            for line in pch:
                line = line[0:72]

                # check for new entity ID
                if line.startswith('$'):
                    entityID = int([i for i in line.split(' ') if i != ''][2])

                    # creates domain and range dictionary for entityID
                    self.data[entityID] = {'domain': [], 'range': []}
                    self.entitylist.append(entityID)

                # add data points to the entity
                elif line.startswith(' '):
                    # parse line
                    datapoint = [float(i) for i in line.split(' ') if i != '']
                    del datapoint[0]

                    #add data to class data attribute
                    self.data[entityID]['domain'].append(datapoint[0])
                    self.data[entityID]['range'].append(datapoint[1])

    def get_entity_list(self):
        """Get a list of each entity (node, element, etc.) in the simple
        punch file."""

        return(self.entitylist)

    def get_domain(self, entityID):
        """Get a list of all values in the domain of the simple punch file"""

        return(self.data[entityID]['domain'])

    def get_range(self, entityID):
        """Get a list of all values in the range of the simple punch file"""

        return(self.data[entityID]['range'])

    def get_plot(self, entityID, fig_num=1, xscale='linear', yscale='linear'):
        """Plot data for single entity using matplotlib.pyplot

        Returns:
            matplotlib.pyplot.plot object

        """

        # create figure
        plt.figure(fig_num)

        # define domain and range
        plot_domain = self.get_domain(entityID)
        plot_range = self.get_range(entityID)

        # plot domain and range to figure
        plt.plot(plot_domain, plot_range, label=entityID)

        # setup plot
        plt.legend(loc='upper left')
        plt.figure(fig_num).add_subplot(1,1,1).set_xscale(xscale)
        plt.figure(fig_num).add_subplot(1,1,1).set_yscale(yscale)

        # show figure
        return plt.figure(fig_num)

    def get_plot_all(self, fig_num, xscale='linear', yscale='linear'):
        """Plot data for all entities on single plot using
        matplotlib.pyplot

        Returns:
            matplotlib.pyplot.plot object

        """

        for each in self.get_entity_list():
            self.get_plot(each, fig_num)

        # setup plot
        plt.legend(loc='upper left')
        plt.figure(fig_num).add_subplot(1,1,1).set_xscale(xscale)
        plt.figure(fig_num).add_subplot(1,1,1).set_yscale(yscale)

    def export(self, targetfile):
        """Export data to CSV file"""

        with open(targetfile, 'w', newline='') as csvfile:

            fields = ['entity', 'domain', 'range']
            writer = csv.DictWriter(csvfile,
                                    fieldnames=fields,
                                    dialect='excel',
                                    delimiter=',')

            writer.writeheader()

            for entity in self.get_entity_list():
                d = self.get_domain(entity)
                r = self.get_range(entity)

                # write in each row
                for i, j in zip(d, r):
                    writer.writerow({'entity': entity,
                                     'domain': i,
                                     'range': j})
