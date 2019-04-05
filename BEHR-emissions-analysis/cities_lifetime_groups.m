classdef cities_lifetime_groups
    %CITIES_LIFETIME_GROUPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant = true)

    end
    
    methods(Static = true)
        function ky = key_years()
            % An empty array = no trend. NaN = all bad
            weekdays.Albuquerque = [2006, 2010];
            weekdays.Jacksonville = [2011, 2012];
            weekdays.Minneapolis = [2006, 2013];
            weekdays.Charlotte = [2008, 2009];
            weekdays.Cincinnati = [2010, 2011];
            weekdays.Indianapolis = [2006, 2013];
            weekdays.Nashville = [2006, 2013];
            
            weekdays.Columbus = [2008, 2009];
            weekdays.Omaha = [2007, 2013];
            weekdays.Pittsburgh = [2011, 2013];
            weekdays.Reno = [2006, 2011];
            weekdays.SanFrancisco = [2006, 2011];
            weekdays.Toronto = [2008, 2013];
            weekdays.Tucson = [2009, 2013];
            weekdays.Orlando = [2007, 2009];
            weekdays.SanDiego = [2006, 2013];
            
            weekdays.Cleveland = [2006, 2009, 2010];
            weekdays.Dallas = [2006, 2010, 2013];
            weekdays.LosAngeles = [2006, 2010, 2013];
            weekdays.Memphis = [2006, 2010, 2013];
            
            weekdays.Knoxville = [2006, 2008, 2012];
            weekdays.Philadelphia = [2006, 2007, 2013];
            weekdays.WashingtonDC = [2006, 2009, 2013];
            
            weekdays.Atlanta = [2006, 2008, 2011];
            weekdays.Denver = [2006, 2010, 2013];
            weekdays.KansasCity = [2006, 2009, 2011];
            weekdays.Montreal = [2006, 2009, 2013];
            weekdays.NewOrleans = [2007, 2009, 2012];
            weekdays.NewYork = [2006, 2011, 2013];
            weekdays.StLouis = [2006, 2007, 2013];
            weekdays.Tampa = [2006, 2011, 2013];
            
            weekdays.Chicago = [2006, 2007, 2010, 2012];
            weekdays.Detroit = [2006, 2009, 2012, 2013];
            weekdays.Richmond = [2006, 2008, 2011, 2013];
            
            weekdays.Houston = [];
            weekdays.SaltLakeCity = [];
            weekdays.Seattle = [];
            
            weekdays.Austin = NaN;
            weekdays.Bakersfield = NaN;
            weekdays.Baltimore = NaN;
            weekdays.Boston = NaN;
            weekdays.Cheyenne = NaN;
            weekdays.Fresno = NaN;
            weekdays.LasVega = NaN;
            weekdays.Miami = NaN;
            weekdays.Phoenix = NaN;
            weekdays.Portland = NaN;
            weekdays.Sacramento = NaN;
            weekdays.SanAntonio = NaN;
            
            
            weekends.Albuquerque = [2006, 2008, 2011, 2012];
            weekends.Atlanta = [2006, 2008, 2013];
            weekends.Austin = NaN;
            weekends.Bakersfield = NaN;
            weekends.Baltimore = NaN;
            weekends.Boston = NaN;
            weekends.Charlotte = [2006, 2009, 2012];
            weekends.Cheyenne = NaN;
            weekends.Chicago = [2006, 2009, 2011, 2013];
            weekends.Cincinnati = [2009 2013];
            weekends.Cleveland = [2006, 2013];
            weekends.Columbus = [2006, 2008, 2010];
            weekends.Dallas = [2006, 2010, 2011, 2012];
            weekends.Denver = [];
            weekends.Detroit = [2006, 2008, 2010, 2013];
            weekends.Fresno = NaN;
            weekends.Houston = [2006, 2007, 2009];
            weekends.Indianapolis = [2006, 2009];
            weekends.Jacksonville = [];
            weekends.KansasCity = [2007, 2010, 2013];
            weekends.Knoxville = [2006, 2011];
            weekends.LasVegas = NaN;
            weekends.LosAngeles = [2006, 2010, 2013];
            weekends.Memphis = [2006, 2008, 2013];
            weekends.Miami = NaN;
            weekends.Minneapolis = [2006, 2008, 2009];
            weekends.Montreal = [2006, 2008, 2012, 2013];
            weekends.Nashville = [2006, 2009, 2013];
            weekends.NewOrleans = [2008, 2010, 2012];
            weekends.NewYork = [2006, 2009, 2013];
            weekends.Omaha = [2007, 2013];
            weekends.Orlando = [2006, 2013];
            weekends.Philadelphia = [2006, 2011, 2013];
            weekends.Phoenix = NaN;
            weekends.Pittsburg = [2006, 2009, 2013];
            weekends.Portland = [2006, 2011, 2013];
            weekends.Reno = [2006, 2011];
            weekends.Richmond = [2006, 2010, 2011];
            weekends.Sacramento = [2007, 2012];
            weekends.SaltLakeCity = [2006, 2011, 2013];
            weekends.SanAntonio = NaN;
            weekends.SanDiego = [2006, 2013];
            weekends.SanFrancisco = NaN;
            weekends.Seattle = [];
            weekends.StLouis = [2006, 2011];
            weekends.Tampa = [2006, 2007, 2010];
            weekends.Toronto = NaN;
            weekends.Tucson = NaN;
            weekends.WashingtonDC = [2006, 2008, 2009, 2012, 2013];
            
            ky = struct('behr_weekdays', weekdays, 'behr_weekends', weekends);
        end
        
        function kyrs = get_key_years(city, days_of_week, varargin)
            % GET_KEY_YEARS( CITY, DAYS_OF_WEEK ) Return an array of the
            % key years that define the lifetime curve for CITY.
            % DAYS_OF_WEEK must be 'TWRF' (weekdays) or 'US' (weekends).
            %
            % GET_KEY_YEARS( CITY, DAYS_OF_WEEK, 'no_nans' ) will return an
            % empty array rather than an NaN if all of the years are
            % invalid fits.
            E = JLLErrors;
            p = advInputParser;
            p.addFlag('no_nans');
            p.parse(varargin{:});
            pout = p.Results;
            
            no_nans = pout.no_nans;
            
            kyrs = cities_lifetime_groups.key_years;
            switch lower(days_of_week)
                case 'twrf'
                    kyrs = kyrs.behr_weekdays;
                case 'us'
                    kyrs = kyrs.behr_weekends;
                    warning('Weekend key years have not been updated since switching to the more conservative t-test')
                otherwise
                    E.badinput('DAYS_OF_WEEK must be "TWRF" or "US"')
            end
            city = regexprep(city, '\s+', '');
            kyrs = kyrs.(city);
            if all(isnan(kyrs)) && no_nans
                kyrs = [];
            end
        end
        
        function cities = decr_lifetime(include_short)
            % DECR_LIFETIME( ) 
            % DECR_LIFETIME( INCLUDE_SHORT ) Return the list of cities with 
            % decreasing lifetime. Set INCLUDE_SHORT to true or omit to
            % include cities with < 60 points in their line densities.
            cities = {'Albuquerque', 'Jacksonville', 'Minneapolis'};
            if nargin < 1 || include_short
                cities = veccat(cities, {'Charlotte', 'Cincinnati', 'Indianapolis', 'Nashville'});
            end
        end
        
        function cities = incr_lifetime(include_short)
            % INCR_LIFETIME( ) 
            % INCR_LIFETIME( INCLUDE_SHORT ) Return the list of cities with 
            % increasing lifetime. Set INCLUDE_SHORT to true or omit to
            % include cities with < 60 points in their line densities.
            cities = {'Pittsburgh', 'Reno', 'San Francisco', 'Toronto', 'Tucson'};
            %cities = {'Pittsburgh', 'Reno', 'Toronto', 'Tucson'};
            if nargin < 1 || include_short 
                cities = veccat(cities, {'Columbus','Orlando', 'San Diego', 'Tampa'});
            end
        end
        
        function cities = ccup_lifetime(include_short)
            % CCUP_LIFETIME( ) 
            % CCUP_LIFETIME( INCLUDE_SHORT ) Return the list of cities with 
            % concave up lifetime. Set INCLUDE_SHORT to true or omit to
            % include cities with < 60 points in their line densities.
            cities = {'Dallas', 'Los Angeles', 'Memphis', 'Omaha'};
            if nargin < 1 || include_short 
                cities = veccat(cities, {'Cleveland','Knoxville', 'Philadelphia', 'Washington DC'});
            end
        end
        
        function cities = ccdown_lifetime(include_short)
            % CCDOWN_LIFETIME( ) 
            % CCDOWN_LIFETIME( INCLUDE_SHORT ) Return the list of cities 
            % with concave down lifetime. Set INCLUDE_SHORT to true or omit 
            % to include cities with < 60 points in their line densities.
            cities = {'Atlanta', 'Denver', 'Kansas City', 'Montreal', 'New Orleans', 'New York', 'St Louis'};
            if nargin < 1 || include_short 
                cities = veccat(cities, {});
            end
        end
        
        function cities = complex_lifetime(include_short)
            % COMPLEX_LIFETIME( ) 
            % COMPLEX_LIFETIME( INCLUDE_SHORT ) Return the list of cities
            % with lifetime changes that do not fit in the other four
            % groups. Set INCLUDE_SHORT to true or omit to include cities
            % with < 60 points in their line densities.
            cities = {'Chicago', 'Detroit'};
            if nargin < 1 || include_short 
                cities = veccat(cities, {'Richmond'});
            end
        end
        
        function cities = no_trend_lifetime(include_short)
            % NO_TREND_LIFETIME( ) 
            % NO_TREND_LIFETIME( INCLUDE_SHORT ) Return the list of cities
            % with no trend in lifetime. Set INCLUDE_SHORT to true or omit
            % to include cities with < 60 points in their line densities.
            cities = {'Houston', 'Salt Lake City'};
            if nargin < 1 || include_short 
                cities = veccat(cities, {'Seattle'});
            end
        end
        
        function cities = all_bad_fit_lifetime(include_short)
            % DECR_LIFETIME( ) Return the list of cities with all bad
            % weekday fits. Set INCLUDE_SHORT to true or omit to include
            % cities with < 60 points in their line densities.
            cities = {'Austin', 'Bakersfield', 'Baltimore', 'Boston', 'Cheyenne', 'Fresno', 'Las Vegas', 'Miami', 'Phoenix', 'Portland', 'Sacramento', 'San Antonio'};
        end
    end
end

