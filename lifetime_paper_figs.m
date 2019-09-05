function lifetime_paper_figs

do_save = false;
do_close = false;

mydir = fileparts(mfilename('fullpath'));
images_dir = fullfile(mydir, '..', '..', 'Science', 'Response', 'Images');
supp_dir = fullfile(images_dir, 'Supplement');

theoretical_lifetime_plot();
ld_example_plot();
lifetime_groups_plot();
no2_vcd_boxplot();
expected_vcds_plot();
lifetime_ens_plot();
weekend_weekday_diff();
make_t_table();

    function theoretical_lifetime_plot
        VOCr = [1, 5, 10];
        fig_main = figure;
        ax_main = gca;
        fig_supp = figure;
        ax_supp = gca;
        fig_species = figure;
        ax_species = gca;
        colors = colormap('lines');
        l_main = gobjects(numel(VOCr),1);
        l_supp = gobjects(numel(VOCr),1);
        
        for i_vocr = 1:numel(VOCr)
            [tau, tau_hno3, tau_ans, nox, species] = nox_lifetime(logspace(-11,-8,70)*2e19, 'vocr', VOCr(i_vocr));
            % plot with NOx as NO2 VCDs assuming a 1 km PBL, a 6 K/km lapse
            % rate, a 4:1 NO2:NO ratio, and a 5e14 free trop background
            nox_vcd = 0.8*(nox .* 1e5) + 5e14;
            l_main(i_vocr) = line(nox_vcd, tau, 'linewidth', 2, 'parent', ax_main, 'color', colors(i_vocr, :));
            
            % plot all three lifetimes for supplement just in case
            l_supp(i_vocr) = line(nox_vcd, tau, 'linestyle', '-', 'linewidth', 2, 'parent', ax_supp, 'color', colors(i_vocr, :));
            line(nox_vcd, tau_hno3, 'linestyle', '--', 'linewidth', 2, 'parent', ax_supp, 'color', colors(i_vocr, :));
            line(nox_vcd, tau_ans, 'linestyle', ':', 'linewidth', 2, 'parent', ax_supp, 'color', colors(i_vocr, :));
            
            % plot the three HOx species
            l_species(1) = line(nox_vcd, species.oh/2e19*1e12, 'linestyle', '-', 'linewidth', 2, 'parent', ax_species, 'color', colors(i_vocr, :));
            l_species(2) = line(nox_vcd, species.ho2/2e19*1e12, 'linestyle', '--', 'linewidth', 2, 'parent', ax_species, 'color', colors(i_vocr, :));
            l_species(3) = line(nox_vcd, species.ro2/2e19*1e12, 'linestyle', ':', 'linewidth', 2, 'parent', ax_species, 'color', colors(i_vocr, :));
        end
        
        voc_lstr = sprintfmulti('VOC_R = %d s^{-1}', VOCr);
        legend(l_main, voc_lstr, 'Location', 'northwest');
        legend(l_supp, voc_lstr, 'Location', 'northwest');
        legend(l_species, {'OH', 'HO_2', 'RO_2'}, 'Location', 'northwest');
        
        all_axes = [ax_main, ax_supp, ax_species];
        ylabel_strs = {'\tau (h)', '\tau (h)', 'Concentration (ppt)'};
        redo_yticks = [false, false, true];
        for a=1:numel(all_axes)
            % Sometimes even plotting with log-log, the x- and y- scale revert
            % to linear.
            set(all_axes(a),'fontsize',14,'xscale','log','yscale','log');
            % It looks nicer to plot the x-ticks as 0.1, 1.0, 10.0, etc.
            % rather than 10^-1, 10^0, etc.
            set(all_axes(a),'XTickLabel', sprintfmulti('%.1f', get(all_axes(a),'XTick')/1e15));
            if redo_yticks(a)
                set(all_axes(a), 'YTickLabel', sprintfmulti('%.1f', get(all_axes(a),'YTick')));
            end
            xlabel(all_axes(a), 'NO_2 VCD (10^{15} molec. cm^{-2})');
            ylabel(all_axes(a), ylabel_strs{a});
        end
        
        % helps ensure a consistent figure size across different machines
        reset_fig_size([fig_main, fig_supp, fig_species]);
        
        if do_save
            save_the_fig(fig_main, 'Theoretical-Lifetime', false);
            save_the_fig(fig_supp, 'HNO3-ANs-Lifetimes', true);
            save_the_fig(fig_species, 'HOx-Species', true);
        end
        if do_close
            close([fig_main, fig_supp, fig_species]);
        end
    end

    function ld_example_plot()
        loc_names = {'Minneapolis', 'Reno', 'St Louis', 'Washington DC'};
        order = [1, 2, 4, 3];
        figs = misc_emissions_analysis.plot_line_dens_fits_by_year('loc_inds', loc_names, 'include_weekends', false,...
            'key_years_only', true, 'ynorm', 'squeeze', 'xnorm', 'fitmax');
        combo_fig = combine_plots(figs(order), 'dims', [2 2], 'scale', 1);
        %close(figs)
        label_subfigs(combo_fig, 'capital', true);
        if do_save
            save_the_fig(combo_fig, 'line-dens-examples', false);
        end
    end

    function no2_vcd_boxplot()
        fig = misc_emissions_analysis.plot_box_city_group_vcds();
        if do_save
            save_the_fig(fig, 'no2-vcds-boxplots-all', false);
        end
    end

    function lifetime_groups_plot()
        fig = figure;
        subplot_stretch(4,1.3);
        
        decr_cities = cities_lifetime_groups.decr_lifetime;
        incr_cities = cities_lifetime_groups.incr_lifetime;
        ccu_cities = cities_lifetime_groups.ccup_lifetime;
        ccd_cities = cities_lifetime_groups.ccdown_lifetime;
        cities = {decr_cities, incr_cities, ccu_cities, ccd_cities};
        
        common_opts = {'plot_quantity', 'Lifetime', 'normalize', true,...
            'always_restrict_to_moves', false, 'req_most', false,...
            'req_num_pts', false, 'min_fits_req', 3, 'recalc_err', false};
            
        for i=1:numel(cities)
            ax = subplot(4,1,i);
            misc_emissions_analysis.plot_avg_lifetime_change('locations',cities{i},...
                'plot_averaging', 'Key years', 'incl_err', true, 'ax', ax, common_opts{:});
            misc_emissions_analysis.plot_avg_lifetime_change('locations',cities{i},...
                'plot_averaging', 'Median', 'incl_err', false, 'ax', ax, common_opts{:});
            ylabel('Norm. Lifetime');
            
            % adding the median adds a "data1" line to the legend. Get rid
            % of it.
            ch = find_named_children(ax);
            for i=1:numel(ch)
                if regcmp(ch(i).DisplayName, 'data')
                    ch(i) = [];
                    break
                end
            end
            legend(ch, {ch.DisplayName}, 'location', 'eastoutside');
        end
        label_subfigs(fig, 'xshift', 0.2);
        
        if do_save
            save_the_fig(fig, 'lifetime-groups', false);
        end
    end

    function lifetime_ens_plot()
        fig = figure;
        subplot_stretch(3,2.6);
        
        decr_cities = cities_lifetime_groups.decr_lifetime;
        incr_cities = cities_lifetime_groups.incr_lifetime;
        ccu_cities = cities_lifetime_groups.ccup_lifetime;
        ccd_cities = cities_lifetime_groups.ccdown_lifetime;
        comp_cities = cities_lifetime_groups.complex_lifetime;
        cities = {decr_cities, incr_cities, ccu_cities, ccd_cities, comp_cities};
        
        common_opts = {'plot_quantity', 'Lifetime', 'normalize', false,...
            'always_restrict_to_moves', false,...
            'req_most', false, 'req_num_pts', false};
        for i=1:numel(cities)
            ax = subplot(3,2,i);
            misc_emissions_analysis.plot_avg_lifetime_change('locations',cities{i},...
                'plot_averaging', 'None', 'incl_err', false, 'ax', ax, common_opts{:});
            ylabel('Lifetime (h)');
        end
        label_subfigs(fig, 'xshift', 0.2);
        
        if do_save
            save_the_fig(fig, 'lifetime-ensembles', true);
        end
    end

    function expected_vcds_plot(exclude_vis)
        if nargin < 1
            exclude_vis = false;
        end
        if exclude_vis
            exclude = cities_lifetime_groups.visually_marginal;
            save_name = 'expected-vcds-excl-visual';
        else
            exclude = [];
            save_name = 'expected-vcds';
        end
        % top panel: individual cities' predicted VCDs (normalized?)
        % bottom panel: normalized average MOVES, VCDS, and expected VCDs
        common_opts = {'normalize', false, 'exclude', exclude,...
            'always_restrict_to_moves', true,... 'req_most', true, 'req_num_pts', true, 'incl_err', false};
            'req_most', true, 'req_num_pts', false, 'incl_err', false};
            
        [~,yrs,moves_emis] = misc_emissions_analysis.plot_avg_lifetime_change(...
            'plot_averaging','None','plot_quantity','MOVES','no_fig',true,common_opts{:});
        [~,~,vcds] = misc_emissions_analysis.plot_avg_lifetime_change(...
            'plot_averaging','None','plot_quantity','VCDs','no_fig',true,common_opts{:});
        
        % The normalization in the above function normalizes each city
        % first, then averages. For the predicted columns, I want to get
        % each city's unnormalized so that large cities have the
        % appropriate effect on the average. For consistency, average these
        % other quantities the same way.
        norm_moves_emis = nanmean(moves_emis,1) ./ nanmean(nanmean(moves_emis,1));
        norm_vcds = nanmean(vcds,1) ./ nanmean(nanmean(vcds,1));
        
        fig = figure;
        ax1 = subplot(2,2,1);
        misc_emissions_analysis.plot_avg_lifetime_change(...
            'plot_averaging','None','plot_quantity','MOVES','ax',ax1,common_opts{:});
        set(ax1, 'yscale', 'log', 'ytick', [5e6, 1e7, 5e7], 'yticklabel', {'5 \times 10^6', '10^7', '5 \times 10^7'});
        ylabel(ax1, 'MOVES emissions (Mg NO_x h^{-1})');
        
        ax2 = subplot(2,2,2);
        misc_emissions_analysis.plot_avg_lifetime_change(...
            'plot_averaging','None','plot_quantity','Lifetime','ax',ax2,common_opts{:});
        ylabel(ax2, 'Lifetime (h)');
        
        ax3 = subplot(2,2,3);
        [~,~,predicted_vcds] = misc_emissions_analysis.plot_avg_lifetime_change(...
            'plot_averaging','None','plot_quantity','Expected VCDs','ax',ax3,common_opts{:});
        ylabel(ax3, 'Predicted NO_x column mass (Mg)')
        
        ax4 = subplot(2,2,4);
        l = gobjects(3,1);
        avg_pred_vcds = nanmean(predicted_vcds,1);
        avg_norm_pred_vcds = avg_pred_vcds ./ mean(avg_pred_vcds);
        l(1)=line(ax4, yrs, avg_norm_pred_vcds, 'color', 'k', 'linewidth', 3);
        l(2)=line(ax4, yrs, norm_vcds, 'color', 'b', 'linestyle', '--', 'linewidth', 3);
        l(3)=line(ax4, yrs, norm_moves_emis, 'color', 'r', 'linestyle', '--', 'linewidth', 3);
        legend(ax4, l, {'Predicted NO_x column','BEHR NO_2 columns','MOVES emissions'}, 'location', 'eastoutside');
        
        subplot_stretch(2,2.6);
        
        label_subfigs(fig, 'xshift', 0.2);
        if do_save
            save_the_fig(fig, save_name, false);
        end
    end


    function weekend_weekday_diff()
        common_opts = {'plot_quantity', 'Weekend - weekday lifetime', 'plot_averaging', 'Median', 'normalize', false,...
            'always_restrict_to_moves', false,...
            'req_most', false, 'req_num_pts', false, 'incl_err', false, 'no_fig', true};
        
        cities = cell(4,1);
        cities_names = {'Decreasing','Increasing','CCU','CCD'};
        cities_styles = struct('color', {'k', 'b', [0 0.5 0], 'r'}, 'linewidth', 3);
        cities{1} = cities_lifetime_groups.decr_lifetime;
        cities{2} = cities_lifetime_groups.incr_lifetime;
        cities{3} = cities_lifetime_groups.ccup_lifetime;
        cities{4} = cities_lifetime_groups.ccdown_lifetime;
        
        fig = figure;
        l = gobjects(numel(cities),1);
        for i_grp = 1:numel(cities)
            [~, x, y] = misc_emissions_analysis.plot_avg_lifetime_change(...
                'locations', cities{i_grp}, common_opts{:});
            l(i_grp) = line(x, y, cities_styles(i_grp));
        end
        ylabel('Weekend - weekday lifetime (h)');
        legend(l, cities_names);
        
        if do_save
            save_the_fig(fig, 'weekend-weekday-lifetime-diffs', true);
        end
    end

    function make_t_table()
        csv_file = fullfile(supp_dir, 'key_years_t_scores.csv');
        misc_emissions_analysis.make_t_score_table(csv_file, 'TWRF')
    end

    function save_the_fig(fig, filename, is_supp)
        
        if is_supp
            save_dir = supp_dir;
        else
            save_dir = images_dir;
        end
        
        fullname = fullfile(save_dir, filename);
        save_fig_paper_formats(fig, fullname);
    end

end
