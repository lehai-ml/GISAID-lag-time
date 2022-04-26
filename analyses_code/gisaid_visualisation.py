import pandas as pd
import numpy as np
#chloropleth map
import folium
from folium import plugins
from branca import colormap as cm
import geopandas as gpd

#calmap
from datetime import datetime
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
from typing import List, Union

class Preprocess_GISAID_files:

    @staticmethod
    def read_gisaid_txt(file: str) -> pd.DataFrame:
        """Read Gisaid output from webscrapping function.
    
        Args:
            file (str): file from webscrapping function.
    
        Returns:
            pd.DataFrame: dataframe in the format of
            Virus | Passage | Accession_ID | Collection_Date | Submission_Date | Length | Host | Location | Originating_Lab | Submitting_Lab
        """
        df = pd.read_csv(file,sep='|',header=None,names=['Virus','Passage','Accession_ID',
                                                          'Collection_Date','Submission_Date',
                                                          'Length','Host','Location','Originating_Lab','Submitting_Lab'],encoding='latin-1')
        
        def get_continent(row: pd.Series) -> str:
            #replace the continent names whitespaces
            location = str.split(row,'/')
            return location[0].strip()
    
        def get_country(row: pd.Series) -> str:
            #replace the country names whitespaces
            location = str.split(row,'/')
            return location[1].strip()
        
        df['Country'] = list(map(get_country,df['Location']))
        df['Continent'] = list(map(get_continent,df['Location']))
        
        return df.loc[:,['Virus','Accession_ID','Collection_Date','Submission_Date','Host','Country','Continent']]
    
    @staticmethod
    def read_tsv_file(file: str) -> pd.DataFrame:
        """Read output downloaded from Gisaid. In tsv format.
    
        Args:
            file (str): Output downloaded from Gisaid
    
        Returns:
            pd.DataFrame: dataframe in the format of
            Virus | Passage | Accession_ID | Collection_Date | Submission_Date | Length | Host | Location | Originating_Lab | Submitting_Lab
        """
        df = pd.read_csv(file,sep='\t')
        df = df.loc[:,['strain','gisaid_epi_isl','date','date_submitted','host','country','region']]
        df.columns = ['Virus','Accession_ID','Collection_Date','Submission_Date','Host','Country','Continent']
        df = df.dropna(subset=['Continent','Collection_Date','Submission_Date','Country'],how='any')
        df['Continent'] = list(map(lambda x: x.replace(' ',''), df['Continent']))
        df['Country'] = list(map(lambda x: x.replace(' ',''), df ['Country']))
        return df
    
    @staticmethod
    def get_log_file(file:str) -> pd.DataFrame:
        """Get log files of total number of samples uploaded.
    
        Args:
            file (str): path to txt log file.
    
        Returns:
            pd.DataFrame: dataframe containing |Date_start | Date_end | total. 
        """
        df = pd.read_csv(file,sep='/',header=None,names=['Date_start','Date_end','Total'])
        df['Continent'] = file.split('\\')[1].replace('_','')
        df = df[df['Total']!=0]
        return df

class Barplot:

    @staticmethod
    def groupby_for_barplot(df: pd.DataFrame,
                            groupby:List[Union[str,pd.core.resample.TimeGrouper]]=['Continent'],
                            aggby:dict = {'Lags':['count','mean']},
                            date_start:str='2021-01-01',
                            date_end:str='2022-04-01'):
        """
        

        Parameters
        ----------
        df : pd.DataFrame
            The metadata data frame.
        groupby : List[Union[str,pd.core.resample.TimeGrouper]], optional
            Groupby parameters. Can be either list of strings, or TimeGrouper(used to resample the dates). The default is ['Continent'].
        aggby : dict, optional
            The column by which operation is performed. The default is {'Lags':['count','mean']}.
        date_start : str, optional
            Select the start date of the aggregation/plot. The default is '2021-01-01'.
        date_end : str, optional
            Select the end date of the aggregation/plot. The default is '2022-04-01'.

        Returns
        -------
        agg_df : TYPE
            Aggregated dataframe. Used for Barplot.

        """
        groupby_keys = []
        for key in groupby:
            try:
                if isinstance(key, pd.core.resample.TimeGrouper):
                    key_name = key.key
                    df[key_name] = pd.to_datetime(df[key_name],format='%Y-%m-%d')
                    df = df[(df[key_name]>=date_start) & (df[key_name]<=date_end)]
                    groupby_keys.append(key_name)
                else:
                    groupby_keys.append(key)
            except AttributeError:
                groupby_keys.append(key)
        agg_df = df.groupby(groupby).agg(aggby)
        if isinstance(agg_df.columns,pd.MultiIndex):
            new_cols = [' '.join(col).strip() for col in agg_df.columns.values]
            agg_df = agg_df.droplevel(0,axis=1)
            agg_df.columns = new_cols
        agg_df = agg_df.reset_index()
        agg_df = agg_df.sort_values(by=groupby_keys)
        agg_df = agg_df.dropna().reset_index(drop=True)
        return agg_df
    
    @staticmethod
    def convert_to_category(data:pd.Series,date_start:str,date_end:str):
        
        return pd.Categorical(data,
                              pd.date_range(date_start,date_end, freq='M').strftime('%Y-%b'),
                              ordered=True)
        
    
    @staticmethod
    def plot(df:pd.DataFrame,
             separateby:str='Continent',
             colorby:str='Lags mean',
             x:str='Collection_Month',
             y:str='Lags count',
             yscalelog:bool=True,
             figsize = (20,10),
             **figkwargs):
        """
        Parameters
        ----------
        df : pd.DataFrame
            See groupby_for_barplot.
        separateby : str, optional
            Separate the subplots by this parameter. The default is 'Continent'.
        colorby : str, optional
            Color the bars by this parameter. The default is 'Lags mean'.
        x : str, optional
            X axis values/ must be categorical data. The default is 'Collection_Month'.
        y : str, optional
            Y axis values/ must be numerical. The default is 'Lags count'.
        yscalelog : bool, optional
            If Yaxis must be scaled by log. The default is True.
        title : str, optional
            Name of the plot title. The default is None.
        **figkwargs:
            xticklabels
            xlable_fontdict
            xlabel
            ylabel
            xlabelpos
            ylabelpos
            title
        Returns
        -------
        None.

        """
        
        total = df[separateby].unique()
        if len(total)>3:
            row = int(np.ceil(len(total)/3))
            column = 3
        else:
            row = 1
            column = len(total)
        fig, axes = plt.subplots(row,column,sharex=True,sharey=True,figsize=figsize)
        all_lags = df[colorby].values
        pal = sns.color_palette('coolwarm', len(all_lags))
        rank = all_lags.argsort().argsort()
        my_cmap = ListedColormap(pal)
        norm = plt.Normalize(all_lags.min(), all_lags.max())
        sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=norm)
        sm.set_array([])
        regions = df[separateby].unique()
        try:
            axes = axes.flatten()
        except AttributeError:
            axes = [axes]
        if yscalelog:
            plt.yscale('log')
        for idx, ax in enumerate(axes):
            data = df[
                df[separateby] == regions[idx]]
            index = data.index.values
            palette = np.array(pal)
            g = sns.barplot(x=x,
                            y=y,
                            data=data,
                            palette=palette[rank[index]],
                            ax=ax)
            ax.set_title(f'{separateby} = {regions[idx]}',fontsize=20)
            ax.set(xlabel=None)
            ax.set(ylabel=None)
            if row > 1:
                if (idx >= 3):
                    if 'xticklabels' not in figkwargs:
                        figkwargs['xticklabels'] = g.get_xticklabels()
                    if 'xlabel_fontdict' not in figkwargs:
                        figkwargs['xlabel_fontdict'] = {'fontsize':10}
                    g.set_xticklabels(figkwargs['xticklabels'], rotation=45,fontdict = figkwargs['xlabel_fontdict'])
            else:
                if 'xticklabels' not in figkwargs:
                    figkwargs['xticklabels'] = g.get_xticklabels()
                if 'xlabel_fontdict' not in figkwargs:
                    figkwargs['xlabel_fontdict'] = {'fontsize':10}
                g.set_xticklabels(figkwargs['xticklabels'], rotation=45,fontdict = figkwargs['xlabel_fontdict'])

        fig.subplots_adjust(right=0.80)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.01, 0.7])
        cbar = g.figure.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Median CST Lags (Days)', size=12)
        if 'xlabel' not in figkwargs:
            figkwargs['xlabel'] = f'{x}'
        if 'xlabel_pos' not in figkwargs:
            figkwargs['xlabel_pos'] = (.48,.01)
        if 'ylabel' not in figkwargs:
            figkwargs['ylabel'] = 'Total sample collected'
        if 'ylabel_pos' not in figkwargs:
            figkwargs['ylabel_pos'] = (.1,.5)
        if 'title' not in figkwargs:
            figkwargs['title'] = None
            
        fig.text(figkwargs['xlabel_pos'][0], figkwargs['xlabel_pos'][1], figkwargs['xlabel'], ha='center', size='xx-large')
        fig.text(figkwargs['ylabel_pos'][0], figkwargs['ylabel_pos'][1],
                 figkwargs['ylabel'],
                 va='center',
                 rotation='vertical',
                 size='xx-large')
        
        fig.suptitle(
            figkwargs['title'],
            size='xx-large')
        
        return fig

class Preprocess_Genbank_files:
    @staticmethod
    def read_csv(self,file):
        return pd.read_csv(file,sep='\t',index_col=[0])
    
    @staticmethod
    def preprocess_genbank_df(self,df):
        """
        Select relevant columns
        
        """
        df = df.loc[:,['strain','gisaidEpiIsl','date','dateSubmitted','host','country','region','pangoLineage']]
        df.columns = ['Virus','Accession_ID','Collection_Date','Submission_Date','Host','Country','Continent','PangoLineage']
        df = df.dropna(subset=['Continent','Collection_Date','Submission_Date','Country'],how='any')
        df['Submission_Date'] =  pd.to_datetime(df['Submission_Date'])
        df['Collection_Date'] =  pd.to_datetime(df['Collection_Date'])
        df['Lags'] = df['Submission_Date'] -  df['Collection_Date']
        df['Collection_Month'] =  df['Collection_Date'].dt.strftime('%Y-%b')
        df['Lags'] = df['Lags'].dt.days
        return df

class Lagmap:
    def __init__(self,all_countries:pd.DataFrame):
        """
        Plot choropleth map: visualise over time across multiple countries

        Parameters
        ----------
        all_countries : pd.DataFrame
            Must have format where the first two columns are groupedby index
                                          Lags
            Country     Collection_Date       
            Afghanistan 2020-05-30       261.0
                        2020-06-02       258.0
        Returns
        -------
        None

        """
        self.all_countries = all_countries.reset_index()
        self.all_countries.columns = ['Country','Date','Value']
    
    def prepare_choropleth_data(self,
                                start:str ='2020-01-01',
                                end:str ='2022-01-01',
                                freq:str ='M',
                                strformat:str='%Y-%b',
                                interpolate:bool =True,
                                plot_all:bool =True):
        """
        

        Parameters
        ----------
        start : str, optional
            DESCRIPTION. The default is '2020-01-01'.
        end : str, optional
            DESCRIPTION. The default is '2022-01-01'.
        freq : str, optional
            DESCRIPTION. The default is 'M'.
        strformat : str, optional
            DESCRIPTION. The default is '%Y-%b'.
        interpolate : bool, optional
            DESCRIPTION. The default is True. If False, make sure the data has no NaN.
        plot_all : bool, optional
            whether to plot all countries. The default is True.
            if False, plot only ones that has missing less than 2 cells.

        Returns
        -------
        None.

        """
        
        self.all_countries['Date'] = pd.Categorical(
            self.all_countries['Date'],
            categories=pd.date_range(start, end, freq=freq).strftime(strformat),
            ordered=True) # order the dates
        
        self.all_countries = self.all_countries.dropna() # remove all the dates that is not in the specified range.
        date_order = self.all_countries['Date'].cat.categories
        all_index = pd.MultiIndex.from_product(
            [self.all_countries['Country'].unique(), date_order],
            names=['Country', 'Date']) # create all possible indices
        self.all_countries = self.all_countries.set_index([
            'Country', 'Date'
        ]).reindex_like(pd.DataFrame({'Value': 0}, index=all_index)) #combine so there is a new dataframe

        self.all_countries = self.all_countries.unstack()  # columns = dates, rows = countries.
        self.all_countries = self.all_countries.droplevel(
            0, axis=1)  # removes the multiindex columns
        self.all_countries = self.all_countries[date_order] # order the columns
        self.all_countries.columns = pd.DatetimeIndex(self.all_countries.columns) # make the columns into a datetimeIndex.
        empty_list = self.get_empty_list(self.all_countries) # get the countries where there is missing date 
        if not plot_all:
            countries_to_plot = [k for k, v in empty_list.items() if v < 2] # v is the number of cell that the corresponding row is missing.
            self.all_countries = self.all_countries.loc[countries_to_plot]
        if interpolate:  # interpolate using available information of time.
            # first_interpolation = self.all_countries.apply(
            #     (lambda x: x.interpolate(method='linear',
            #                               limit_direction='both')),
            #     axis=0)
            self.all_countries = self.all_countries.apply(
                (lambda x: x.interpolate(method='time', limit_direction='both')
                  ),
                axis=1)        

        self.all_countries = self.all_countries.stack(dropna=False)
        self.all_countries = pd.DataFrame(self.all_countries.rename('Value'))

    def get_style_dict(self, all_countries: pd.DataFrame, geodata:dict):
        """
        

        Parameters
        ----------
        all_countries : pd.DataFrame
            Must have format where the first two columns are groupedby index
                                          Lags
            Country     Collection_Date       
            Afghanistan 2020-05-30       261.0
                        2020-06-02       258.0
        geodata : dict
            json loaded data containing the data of the geographical locations.

        Returns
        -------
        countries_gdf : gpd.GeoDataFrame
            DESCRIPTION.
        style_dict : dict
            DESCRIPTION.
        cmap : list
            DESCRIPTION.

        """
        
        all_countries = all_countries.reset_index()
        all_countries['Date'] = pd.DatetimeIndex(
            all_countries['Date'])
        all_countries = rename_countries(all_countries)

        # convert it to unix nanoseconds
        all_countries['Date_sec'] = all_countries['Date'].view(
            np.int64) // 10**9
        all_countries['Date_sec'] = all_countries['Date_sec'].astype(str)

        geodata = gpd.GeoDataFrame.from_features(geodata)
        geodata.columns = ['geometry', 'Country']
        all_countries = pd.merge(all_countries, geodata, on='Country')

        max_colour = max(all_countries['Value'])
        min_colour = min(all_countries['Value'])
        cmap = cm.linear.YlOrBr_09.scale(min_colour, max_colour)
        all_countries['colour'] = all_countries['Value'].map(cmap)
        country_list = all_countries['Country'].unique().tolist()
        country_idx = range(len(country_list))

        style_dict = {}
        for i in country_idx:
            country = country_list[i]
            result = all_countries[all_countries['Country'] == country]
            inner_dict = {}
            for _, r in result.iterrows():
                inner_dict[r['Date_sec']] = {
                    'color': r['colour'],
                    'opacity': 0.7
                }
            style_dict[str(i)] = inner_dict

        countries_gdf = gpd.GeoDataFrame(all_countries[[
            'geometry'
        ]]).drop_duplicates().reset_index(drop=True)

        return countries_gdf, style_dict, cmap

    def plot_choropleth(self, geodata:dict):
        """
        

        Parameters
        ----------
        geodata : dict
            json loaded coordinates of the countries.

        Returns
        -------
        m : folium.Map
            choropleth map.

        """
        

        countries_gdf, style_dict, cmap = self.get_style_dict(self.all_countries,
                                                              geodata=geodata)
        m = folium.Map([50.0, 5.0],
                       tiles="openstreetmap",
                       zoom_start=4,
                       max_bounds=True)
        _ = plugins.TimeSliderChoropleth(data=countries_gdf.to_json(),
                                         styledict=style_dict).add_to(m)

        _ = cmap.add_to(m)
        return m

    def get_empty_list(self, df):
        """
        Get the countries where the dates are missing
        """
        empty_list = {}
        try:
            df_unstacked = df.droplevel(0, axis=1)
        except ValueError:
            df_unstacked = df
        for country in df_unstacked.index:
            empty_list[country] = len(
                df_unstacked.columns[(df.loc[country, :].isnull())].to_list())

        return empty_list

def rename_countries(df):
    df = df.replace('UnitedKingdom', 'United Kingdom')
    df = df.replace('CzechRepublic', 'Czech Republic')
    df = df.replace('BosniaandHerzegovina', 'Bosnia and Herzegovina')
    df = df.replace('BurkinaFaso', 'Burkina Faso')
    df = df.replace('DemocraticRepublicoftheCongo',
                    'Democratic Republic of the Congo')
    df = df.replace('Congo (Democratic Republic)',
                    'Democratic Republic of the Congo')
    df = df.replace('Congo',
                    'Republic of the Congo')
    df = df.replace('CaboVerde','Cabo Verde')
    df = df.replace('Cape Verde','Cabo Verde')
    df = df.replace('Czechia','Czech Republic')
    df = df.replace("Cote d'Ivoire",'Ivory Coast')
    df = df.replace('Curaçao','Curacao')
    df = df.replace('CostaRica','Costa Rica')
    df = df.replace('DominicanRepublic', 'Dominican Republic')
    df = df.replace('Dominica', 'Dominican Republic')
    df = df.replace('ElSalvador', 'El Salvador')
    df = df.replace('EquatorialGuinea', 'Equatorial Guinea')
    df = df.replace('Guinea-Bissau', 'Guinea Bissau')
    df = df.replace('Gambia, The', 'Gambia')
    df = df.replace('HongKong', 'Hong Kong')
    df = df.replace('Myanmar (Burma)', 'Myanmar')
    df = df.replace('NewZealand', 'New Zealand')
    df = df.replace('NorthAmerica','North America')
    df = df.replace('NorthMacedonia', 'Macedonia')
    df = df.replace('PapuaNewGuinea', 'Papua New Guinea')
    df = df.replace('PuertoRico','Puerto Rico')
    df = df.replace('RepublicoftheCongo', 'Republic of the Congo')
    df = df.replace('SaudiArabia', 'Saudi Arabia')
    df = df.replace('SierraLeone', 'Sierra Leone')
    df = df.replace('SouthAfrica', 'South Africa')
    df = df.replace('SouthKorea', 'South Korea')
    df = df.replace('Korea (South)', 'South Korea')

    df = df.replace('SouthSudan', 'South Sudan')
    df = df.replace('SriLanka', 'Sri Lanka')
    df = df.replace('Réunion',
                    'Reunion')
    df = df.replace('TheBahamas', 'The Bahamas')
    df = df.replace('Bahamas, The', 'The Bahamas')
    df = df.replace('Timor-Leste', 'East Timor')
    df = df.replace('TrinidadandTobago', 'Trinidad and Tobago')
    df = df.replace('USA', 'United States of America')
    df = df.replace('United States', 'United States of America')
    df = df.replace('US Virgin Islands', 'U.S. Virgin Islands')
    df = df.replace('Wallis and Futuna', 'Wallis and Futuna Islands')
    df = df.replace('St Lucia',
                    'Saint Lucia')
    df = df.replace('SaintMartin',
                    'St Martin (French part)')
    df = df.replace('St Vincent',
                    'Saint Vincent and the Grenadines')
    
    return df
    
class Calmap:
    def __init__(self,year:int,data:pd.DataFrame):
        """
        

        Parameters
        ----------
        year : int
            DESCRIPTION.
        data : pd.DataFrame
            in the form of: |Date|Value| where date is the index.
            

        Returns
        -------
        data: numpy array
        in the shape 53x7. dates with no value is filled with 0

        """
        
        self.data_array = self.prepare_data_for_calmap(year, data)
        self.year = year
    
    def calmap(self,ax,**kwargs):
        """
        Draw Calendar Heatmap, with borders for each month
        ax = subplot
        """
        ax.tick_params('x',length=0,labelsize='medium',which='major')
        ax.tick_params('y',length=0,labelsize='x-small',which='major')
        
        #creating month borders
        xticks, labels = [],[]
        
        start = datetime(self.year,1,1).weekday()
        
        for month in range(1,13):
            first = datetime(self.year,month,1)
            last = first + relativedelta(months=1,days=-1) #last day for that month
            
            y0 = first.weekday()
            y1 = last.weekday()
            x0 = (int(first.strftime('%j')) + start -1 )//7
            x1 = (int(last.strftime('%j')) + start -1 )//7
            
            #define border
            P = [(x0,y0),
                (x0,7),
                (x1,7),
                (x1,y1+1),
                (x1+1,y1+1),
                (x1+1,0),
                (x0+1,0),
                (x0+1,y0)]
            
            xticks.append(x0+(x1-x0+1)/2)
            labels.append(first.strftime('%b'))
            poly = Polygon(P,edgecolor = 'black',facecolor ='None', linewidth=1,zorder=20,clip_on=False)
            ax.add_artist(poly)
            
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(0.5+ np.arange(7))
        ax.set_yticklabels(['Mon','Tue','Wed','Thu','Fri','Sat','Sun'])
        ax.set_title('{}'.format(self.year),weight='semibold')
        
        im = ax.imshow(self.data_array,extent=[0,53,0,7],zorder=10,origin='lower',cmap='jet',alpha=.75,**kwargs)
        
        return im


    def prepare_data_for_calmap(self,year:int,data:pd.DataFrame)->np.array:
        """
        

        Parameters
        ----------
        year : int
            DESCRIPTION.
        data : pd.DataFrame
            in the form of: |Date|Value| where date is the index.
            

        Returns
        -------
        data: numpy array
        in the shape 53x7. dates with no value is filled with 0

        """
        first = datetime(year,1,1)
        last = datetime(year,12,31)
        all_dates = pd.Series(pd.date_range(first,last).strftime('%Y-%m-%d'),name='Date')
        if isinstance(data.index,pd.DatetimeIndex):
            data.index = data.index.strftime('%Y-%m-%d')
        new_data = data.reindex(all_dates)
        new_data = new_data.fillna(0)
        new_data = new_data.values.flatten()
        
        #add in paddings to make it 53*7
        
        new_data = np.append(np.append(np.zeros((first.weekday(),)),new_data),
                             np.zeros((6-last.weekday(),)))
        new_data[:first.weekday()] = np.nan
        if last.weekday()!=6:
            new_data[(last.weekday()-6):] = np.nan
        new_data = new_data.reshape(53,7)
    #     print(first.weekday(),last.weekday())
        return new_data.T
    
    @staticmethod
    def compare_multiple_maps(axes,year_list,country_list,data,**kwargs):
        for year, ax in zip(year_list,axes):
            for country, ax_column in zip(country_list,ax):
                temp_df = data[(data['Country']==country)
                             & (data['Collection date']<=str(year+1))
                             & (data['Collection date']>str(year))]
                im = Calmap(year, temp_df.set_index('Collection date')[['Lags count']]).calmap(ax_column,**kwargs)
                
                
                
                
            
        
        