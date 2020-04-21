import pandas as pd
import math as mt


class loadcable:
    """Class to loading cable files from a .xslx database

    Returns:
        Dataframe -- Dataframe of the selected cable from database
    """

    def __init__(self, database):
        self.database = database

    def _selectdb(self):
        return pd.read_excel(self.database)

    def filtercable(self, cable_type, cable_name='', cable_awg=None):
        bdcables = self._selectdb()
        cond1 = bdcables.Type == cable_type
        cond2 = bdcables.Name == cable_name

        if cable_awg is None:
            cable = bdcables[cond1 & cond2]
        else:
            cond3 = bdcables.AWG == cable_awg
            cable = bdcables[cond1 & (cond2 | cond3)]
        # pd.set_option('max_colwidth', 33)
        # pd.set_option('expand_frame_repr', False)
        return pd.DataFrame(cable)

    def optimalcable(self, cable_type, cable_name):
        """
        Filter the conductors in a database refereed by the electrical resistance DC

        Returns alternative cables for the original cable chose
        :param cable_type: type of original cable
        :param cable_name: name of original cable
        :return: alternative cables  according electrical resistance 20 of original cable
        """

        cable = {'Type': cable_type, 'Name': cable_name}
        df_cable = self.filtercable(
            cable['Type'], cable['Name'])
        cable_res = df_cable['ElectRes_CC_20'].iloc[0]

        # Show different cables types in database
        bdcables = self._selectdb()
        cables_type_iter = set(bdcables['Type'])
        # Remove the original cable type
        cables_type_iter.remove(cable['Type'])

        # c1 = df_cable?
        # c1 = original cable data row from database
        c1 = bdcables[(bdcables['Type'] == cable['Type']) &
                      (bdcables['Name'] == cable['Name'])]

        # c2 serveral cable types with ElectRes_CC_20 <= original resistance
        # c2 gets only the first cable of each type with ElectRes_CC_20 <= original resistance
        # [:1], and concatenate the rows for each end loop
        c2 = pd.concat([bdcables[(bdcables['Type'] == i) &
                                 (bdcables['ElectRes_CC_20'] <= cable_res)][:1]
                        for i in cables_type_iter], ignore_index=True)

        # optional_cables_df = DataFrame of cables choices for each type, including the original one
        optional_cables_df = pd.concat([c1, c2], ignore_index=True)

        return optional_cables_df
