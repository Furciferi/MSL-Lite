import os

import yaml
import pickle

import pandas as pd
import numpy as np
import math

from astropy import units as u
from astropy.coordinates import SkyCoord,Angle
from astropy.wcs import WCS
from astropy.io import fits

from tqdm import tqdm

from msl_exceptions import MSLSequenceError

print(pd.__version__)
class MSL:
    def __init__(self, debug: bool = False, debug_size: int = 0,
                 low_mem: bool = False):
        """ Initialisation of Master Source List and associated methods"""

        self._config_filepath = "./MSL_config.yaml"
        self._cwd = "./"
        self._master_df = None
        self._args_loaded = False
        self._debug = debug
        self._debug_size = debug_size

        # If needed, clean memory as much as possible if needed
        self._low_mem = low_mem

        self._load_config_file()

        if self._debug:
            tqdm.pandas(ascii=True)


    def _load_config_file(self):
        """
            Method to read and load key value from the parameter file.
        """
        with open(self._config_filepath,'r') as fin:
            try:
                self._params = yaml.load(fin, Loader = yaml.CLoader)
            except AttributeError:
                print("Couldnt use Cloader, using Loader")
                self._params = yaml.load(fin, Loader = yaml.Loader)

    @staticmethod
    def apply_xcs_name(x):
        RA = float(x["RA"])
        DEC = float(x["DEC"])
        RA_h,RA_m,RA_s = Angle(RA*u.degree).hms
        _,DEC_d,DEC_m,DEC_s = Angle(DEC*u.degree).signed_dms
        if DEC>=0:
            sign='+'
        else:
            sign ='-'
        DEC_d=abs(DEC_d)
        xcs_name = ("XMMXCS J{h:02.0f}{r_m:02.0f}{r_s:04.1f}{sign}"
                    "{d:02.0f}{d_m:02.0f}{d_s:04.1f}"
                    ).format(
                                h=RA_h,r_m=RA_m,r_s=RA_s,
                                sign=sign,d=DEC_d,d_m=DEC_m,d_s=DEC_s)
        return xcs_name


    def set_config_filepath(self,path: str):
        """
            Setter for config file path

            :param path: System filepath for the config file
            :type path: str
            :return: None
        """
        if os.path.exists(path):
            self._config_filepath = path
        else:
            print("Config File does not exists: {}".format(path))


    def reload_params(self):
        """
            Reload the params from the config file
        """
        with open(self._config_filepath,'r') as fin:
            try:
                self._params = yaml.load(fin, Loader = yaml.CLoader)
            except AttributeError:
                print("Couldnt use Cloader, using Loader")
                self._params = yaml.load(fin, Loader = yaml.Loader)

    def process_idresults_files(self):
        id_results_dir = self._params["id_results_loc"]
        obsids = next(os.walk(id_results_dir))[1]
        for itter,filename in enumerate(
                        self._params["results_files_to_process"]["in_files"]):
            outfile = self._params["results_files_to_process"]["out_files"][itter]
            print("Processing filenames:",filename, "to",outfile)
            with open(outfile,'w+') as fout:
                for ID in tqdm(obsids,ascii=True):
                    try:
                        int(ID)
                    except:
                        print("Bad OBSID DIR:",ID)
                        continue
                    file_path = os.path.join(   id_results_dir,
                                                ID,
                                                filename)
                    if os.path.isfile(file_path):
                        with open(file_path,'r') as fin:
                            try:
                                next(fin)
                            except:
                                print("Empty file:",file_path)
                                continue
                            for line in fin:
                                print(" ".join([ID,line.rstrip()]),file=fout)


    def aggregate_finalsigs(self):
        """
            Method to aggreagate all finalsig files from each observation.

        """
        if self._debug:
            print("Aggregating Finalsig files")

        det_params = self._params["detections"]
        # Load the finalsigs data in
        sig_file = det_params["sig_file"]
        finalsig_header = det_params["finalsig_header"]
        df_sig = pd.read_csv(   sig_file,sep='\s+',header=None,
                                error_bad_lines=False,
                                warn_bad_lines=True)

        # the file does not have a header so load the param header in
        df_sig.columns = finalsig_header

        self._df_finalsigs = df_sig


    def aggregate_srcsums(self):
        """
            Method to aggregate all srcsums files from each observation.
        """
        if self._debug:
            print("Aggregating Srcsum files")
        det_params = self._params["detections"]

        # Load in all data from the all_sources file
        source_file = det_params["src_file"]
        source_header = det_params["in_header"]
        df_src = pd.read_csv(source_file,sep='\s+',header=None,
                                error_bad_lines=False,
                                warn_bad_lines=True)

        # the file does not have a header so load the param header in
        df_src.columns = source_header

        self._df_srcsums = df_src


    @staticmethod
    def fix_OBSID(x):
        return str(int(x)).zfill(10)

    def create_master_df(self):
        """
            Method to merge finalsigs and srcsums into one master database
            and remove duplicate columns, adding a unique identifier for each
            source.
        """
        if self._debug:
            print("Creating master dataframe from finalsigs and srcsums")
        try:
            df_sigs = self._df_finalsigs
            df_src = self._df_srcsums
        except NameError:
            print("Either finalsigs or srcsums dataframe has not been created")
            raise MSLSequenceError

        def merge_obsid_row(x):
            return x["OBSID"]+"_"+str(x["ROW"])

        # Make sure obsids are 10 digits
        df_sigs["OBSID"] = df_sigs.OBSID.progress_apply(MSL.fix_OBSID)
        df_src["OBSID"] = df_src.OBSID.progress_apply(MSL.fix_OBSID)

        # And provide the unique columns to merge on from obsid and row_num
        df_sigs["UNIQUE"] = df_sigs.progress_apply(merge_obsid_row,axis=1)
        df_src["UNIQUE"] = df_src.progress_apply(merge_obsid_row,axis=1)

        # merge the two files on the unique column, adding "_y" to dupe columns
        df_master = df_src.merge(df_sigs,on="UNIQUE",
                                    suffixes=["","_spare"],how="inner")

        # throw away the dupelicate columns
        to_drop = [x for x in df_master if x.endswith('_spare')]
        df_master.drop(to_drop,axis=1,inplace=True)
        print(df_master.columns)
        # Create a new column containing the XCS naming convention
        df_master.insert(   1,
                            "XCS_NAME",
                            df_master.progress_apply(
                                        MSL.apply_xcs_name,
                                        axis=1),
                            False)

        # Clean up older dataframes if low memory needed
        if self._low_mem:
            del self._df_srcsums
            del self._df_finalsigs
        if self._debug_size:
                df_master = df_master.sample(self._debug_size
                                        ).reset_index(drop=True)
        self._df_master = df_master


    def set_column_float(self, column: str):
        """
            Change datatype of a column to float in the master dataframe
            :param column: Columns in the dataframe to change format
            :type column: Str
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError
        if self._debug:
            print("Setting column ({}) as float".format(column))
        df_master[column].astype(float, copy=False)


    def sort_master_by(self, column: str, ascending: bool = False,
                        return_df: bool = False):
        """
            Method to sort the master dataframe by a column

            :param column: Columns name in which to order the master dataframe
            :type column: Str
            :param ascending: Sort the dataframe ascending, default False
            :type ascending: bool
            :param return_df: Sets whether the master table is returned by the
             function. Default False.
             :type return_df: bool

        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError
        if self._debug:
            print("Sorting Master on {}".format(column))
        df_master.sort_values(  by = [column],
                                axis = 0,
                                inplace = True,
                                ascending = ascending,
                                na_position = "last")

        if return_df:
            return df_master


    def set_master_ID(self):
        """
            Add interger row number as unique identifier to master dataframe
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError
        if self._debug:
            print("Adding master ID to master datagframe")
        # Add ID column to first in master df and not allow if ID already exists
        df_master.insert(0,"ID",range(df_master.shape[0]),False)


    def write_master_to_csv(self, path: str = "",
                            header: bool = True, index: bool = False):
        """
            Write the master dataframe to a csv file
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError
        if path!="":
            df_master.to_csv(path, header = header, index = index)
        else:
            df_master.to_csv(self._params["master_df_csv"], header = header,
                                index = index)


    def load_master_from_csv(self, path: str, index_col = None,
                                fix_obs: bool = True):
        """
            Set master dataframe from a csv file
        """
        if self._debug:
            print("Loading MSL from csv {}".format(path))
        if os.path.exists(path):
            self._df_master = pd.read_csv(path, index_col = index_col)
        else:
            print("CSV file could not be found: {}".format(path))

        self.all_columns_if_float()

        if fix_obs:
            if self._debug:
                print("fixing OBSID pandas float")
            self._df_master["OBSID"] = self._df_master["OBSID"].apply(
                                                    MSL.fix_OBSID)


    @staticmethod
    def get_unique_and_duplicates(dataframe: pd.DataFrame, radius: float,
                                    verbose: bool = False):
        """
            Calculate duplicate sources within a catalog, using astropy's
            SkyCoord.

            :param dataframe: Input catalog with columns RA,DEC in degrees
            :type dataframe: pd.DataFrame
            :param radius: radius to set objects as duplicates, in arcseconds
            :type radius: float
            :return: Unique dataframe, duplicate dataframe
        """
        if verbose:
            print("Converting Dataframe to Astropy Skycoords")
        # Create a catalogue of the sky positions using astropy
        cat = SkyCoord( ra = dataframe.RA.to_numpy()*u.degree,
                        dec = dataframe.DEC.to_numpy() *u.degree)

        if verbose:
            print("Calculating matches against itself - this may take a while")
        # Match the catalog against itself to find obsids
        idxc, idxcatalog, d2d, d3d = cat.search_around_sky(cat,
                                                    radius*u.arcsecond)

        # Work through the matched rows, if the objects are in the same
        # observation, then not a duplicate, but if in another, count as
        # duplicate. Mark IDs already matched, so we can skip them as while
        # assume that the dataframe is ordered with the best sources at the top
        # and the "worse" sources descending.

        seen = np.zeros(dataframe.shape[0])
        dupes = {}
        unique = np.zeros(dataframe.shape[0])
        for index in tqdm(list(set(idxc)),ascii=True):
            if seen[index]:
                # Already said this source is a duplicate
                continue
            seen[index]=1
            # Grab all the idxs which where within the radius
            match_index = idxcatalog[idxc==index]
            # If it matched to itself, no other sources so unique
            if match_index.shape[0]==1:
                unique[index]=1
            else:
                # check if the obsids are different - ie same source
                # if the same obsids, then they are different sources.
                no_clash =[]
                for ind in match_index:
                    if dataframe.loc[index,"OBSID"]!=dataframe.loc[ind,"OBSID"]:
                        no_clash += [dataframe.loc[ind,"ID"]]
                        seen[ind]=1
                # write list of duplicates associated with the "best" source
                # to the duplicate dataframe incase we need to look at these
                # at a later date/analysis
                dupes[dataframe.loc[index,"ID"]]=no_clash
        return dataframe[unique.astype(bool)]["ID"],dupes

    @staticmethod
    def get_source_type(extent: float):
        """
            Method to calculate if a source is extended or point like
        """
        if extent !=0:
            return "Point"
        else:
            return "Ext"
    def flag_source_type(self):
        """
            Add a column to the master dataframe weather the source is a point
            source or extended sources
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        df_master["EXT"] = df_master["PPOINT"].progress_apply(MSL.get_source_type)

    def remove_duplicates(self):
        """
            Reduce master table to unique sources for extended and point like.
            This assumes sources within a given radius are duplicates if they
            appear in other OBSIDS.
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        # Split master into point and extended dataframes
        df_extended = df_master[df_master["EXT"] == 'Ext'].reset_index(drop=True)
        df_point = df_master[df_master["EXT"] == "Point"].reset_index(drop=True)

        if self._debug:
            print("Performing self match on extended - size {}".format(
                            df_extended.shape[0]))
        ext_ids, ext_dupes = MSL.get_unique_and_duplicates(
                                dataframe = df_extended,
                                radius = 30,
                                verbose = self._debug
                                )
        if self._debug:
            print("Performing self match on point sources - size {}".format(
                            df_point.shape[0]))
        pnt_ids, pnt_dupes = MSL.get_unique_and_duplicates(
                                dataframe = df_point,
                                radius = 5,
                                verbose = self._debug
                                )

        # Now we have the list of unique and duplicate IDs we reduce the
        # master dataframe to only the unique IDs
        df_master = df_master[df_master.index.isin(ext_ids + pnt_ids)]

        # TODO:
        # And we create a new dictionary of the duplicate IDs and write to
        # a pickle file.
        # self._duplicates = {**ext_dupes,**pnt_dupes}
        # print(self._duplicates)
        # with open("Duplicates.pkl","wb") as fout:
        #     pickle.dump(self._duplicates, fout, pickle.HIGHEST_PROTOCOL)

    def reformat_flag(self, column: str):
        """
            Replace "ok" with n and otherwise y
        """

        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        def swap(x):
            if x == "ok":
                return "n"
            else:
                return "y"

        df_master[column] = df_master[column].progress_apply(swap)


    def add_OBSID_wcs(self):
        """
            For each unique OBSID, fetch the WCS data from the fits files
            and add 4 columns to the master dataframe, (RA/DEC/PA)_PNT,CDELT2.
        """

        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        # Keys we need for a wcs, add empty columns to master df
        keys = self._params["merged_header"]
        for key in keys:
            df_master[key]=""

        # For each obsid, grab the details and if theyre missing replace the
        # column with a NONE. Insert the value for all rows with that obsid
        obsid_set = list(set(df_master["OBSID"]))
        for obsid in tqdm(obsid_set):
            obsid_path = os.path.join(
                                self._params["data_path"],
                                self._params["fits_path"].format(obs=obsid))
            df_master[df_master["OBSID"]==obsid][key]
            if os.path.exists(obsid_path):
                header = fits.open(obsid_path)[0].header
                for key in keys:
                    try:
                        df_master.loc[df_master["OBSID"]==obsid, key] = str(header[key])
                    except (KeyError, ValueError):
                        df_master.loc[df_master["OBSID"]==obsid, key] = "NONE"

            else:
                print("Obsid Path does not exist:{}".format(obsid_path))
                for key in keys:
                    df_master.loc[df_master["OBSID"]==obsid, key] = "NONE"

    def add_theta_phi_to_master(self):
        """
            Replaces xy2polar in the IDL version. Calcualtes theta and phi
            from OBSID wcs.
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        keys = self._params["merged_header"]

        for key in keys:
            if key not in df_master.columns:
                print(( "Please make sure wcs information is"
                        " contained within the master "
                        "dataframe. Missing"
                        "{}".format(key)))
                raise MSLSequenceError
        list_theta = []
        list_phi = []
        bad_idx = []
        for idx in tqdm(df_master.index):
            header = {}

            # TODO: Remove when OBSIDS are fixed
            if "NONE" in df_master.loc[idx,keys].values:
                bad_idx.append(idx)
                list_theta.append("NONE")
                list_phi.append("NONE")
                continue

            for key in keys:
                header[key] = df_master.loc[idx,key]
            w = WCS(header)

            x_pnt,y_pnt = w.wcs_world2pix(  df_master.loc[idx,"RA_PNT"],
                                            df_master.loc[idx,"DEC_PNT"],1)
            x = df_master.loc[idx,"CENTRAL_X"]
            y = df_master.loc[idx,"CENTRAL_Y"]

            off = abs(y - y_pnt) / abs(x - x_pnt)
            phi = 90 - (180/np.pi)*math.atan(off)
            phi += 360 - df_master.loc[idx, "PA_PNT"]
            phi %= 360

            list_phi.append(phi)

            dist_in_pix = np.sqrt( ( x - x_pnt)**2 + (y - y_pnt)**2)

            list_theta.append(dist_in_pix * df_master.loc[idx,"CDELT2"]*60)

        df_master["PHI"] = list_phi
        df_master["THETA"] = list_theta
        df_master = df_master.drop(bad_idx).reset_index(drop=True)


    def add_sig_error(self):
        """
            Method which uses the PHI value to add fixed sig error values to
            the master dataframe.
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        if "PHI" not in df_master.columns:
            print(("PHI not in master dataframe, please run "
                    "add_theta_phi_to_master first"))
            raise MSLSequenceError
        if self._debug:
            print("Adding Sig Error")

        # load the static values for error from the params file
        sig_x = self._params["detections"]["sig_err_x"]
        sig_y = self._params["detections"]["sig_err_y"]


        xerr = []
        yerr = []

        def calc_sig_err(row,sig_x, sig_y):
            theta = row["THETA"]
            if theta == "NONE":
                return "NONE","NONE"
            elif theta > 15:
                theta = 15
            sig_idx = int(theta)

            return (sig_x[sig_idx], sig_y[sig_idx])

        tmp = df_master.progress_apply(calc_sig_err,args=(sig_x,sig_y), axis=1, result_type='expand')
        df_master["XERR"] = tmp[0]
        df_master["YERR"] = tmp[1]


    def log_sig(self):
        """
            Method which applies -log(sig) to the read in sig value, unless
            0 in which the default is set to -1000.0
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        def logit(sig):
            if sig!=0:
                return -1* math.log(sig)
            else:
                return -1000.0
        df_master.loc[:,"SIG_LOG"] = df_master["SIG"].apply(logit)


    def make_source_csv(self):
        """
            Method to create a csv output file, generally refered to as the
            master source list.
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        with open(self._params["output_csvs"]["source"],'r') as fin:
            try:
                src_params = yaml.load(fin, Loader = yaml.CLoader)
            except AttributeError:
                print("Couldnt use Cloader, using Loader")
                src_params = yaml.load(fin, Loader = yaml.Loader)

        # Add pointless Columns
        df_master["HARD_COUNTS"] = "na"
        df_master["HARD_SIG"] = "na"
        df_master["HARD_RATE"] = "na"
        df_master["HARD_RATIO"] = "err"

        # Add pointless flags
        def cls_flag(x):
            if x == "Ext":
                return "n"
            else:
                return 'e'

        df_master["TARGET_FLAG"] = df_master["EXT"].apply(cls_flag)
        df_master["CLS_FLAG"] = df_master["EXT"].apply(cls_flag)
        df_master["SNR_FLAG"] = df_master["EXT"].apply(cls_flag)

        # And write the MSL
        df_master[src_params["header"]].to_csv(src_params["outfile"],
                                                index=None,header=None)


    def all_columns_if_float(self):
        """
            Simple method to fix pandas reading strings not floats
        """
        try:
            df_master  = self._df_master
        except NameError:
            print("Master dataframe does not exist")
            raise MSLSequenceError

        def try_float(x):
            try:
                x = float(x)
            except ValueError:
                pass
            return x

        for column in tqdm(df_master.columns):
            df_master[column] = df_master[column].apply(try_float)





if __name__ == "__main__":

    msl = MSL(debug=True)
    msl.process_idresults_files()
    msl.aggregate_srcsums()
    msl.aggregate_finalsigs()
    msl.create_master_df()
    msl.set_column_float(column = "SOURCE_COUNTS")
    msl.sort_master_by(column = "SOURCE_COUNTS")
    msl.set_master_ID()
    msl.flag_source_type()
    msl.remove_duplicates()
    msl.sort_master_by(column = "OBSID")
    msl.log_sig()
    msl.make_source_csv()
    msl.write_master_to_csv("./MSL.csv")
