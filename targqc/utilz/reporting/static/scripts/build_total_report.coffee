sampleReport =
    sample:
        name: ''
        display_name: ''
        phenotype: ''
        bam: ''
        bed: ''
        vcf_by_caller:
            name: ''
            summary_qc_rep_fpaths: []
            anno_vcf_fpaths: {}
            anno_filt_vcf_fpaths: {}
    html_fpath: ''
    link: ''
    records: []

record =
    metric: null
    value: ''
    meta: null
    html_fpath: ''

metric =
    name: ''
    short_name: ''
    description: ''
    quality: ''
    presision: 0
    type: null
    all_values_equal: false


DRAGGABLE_COLUMNS = false


################
# Color heat map
BLUE_HUE = 240
BLUE_OUTER_BRT = 55
BLUE_INNER_BRT = 65

GREEN_HUE = 120
GREEN_OUTER_BRT = 50
GREEN_INNER_BRT = 60

RED_HUE = 0
RED_OUTER_BRT = 50
RED_INNER_BRT = 60

MIN_NORMAL_BRT = 80
MEDIAN_BRT = 100  # just white.

get_color = (hue, lightness) ->
    lightness = if lightness? then lightness else 92
    # lightness = Math.round (Math.pow hue - 75, 2) / 350 + 35
    [r, g, b] = hslToRgb(hue / 360, 0.8, lightness / 100)
    return '#' + r.toString(16) + g.toString(16) + b.toString(16)
#    return 'hsl(' + hue + ', 80%, ' + lightness + '%)'

################


check_all_values_equal = (vals) ->
    first_val = null
    for val in vals
        if first_val?
            if val != first_val
                return false
        else
            first_val = val
    return true


get_meta_tag_contents = (rec) ->
    meta = rec.meta

    if meta? and (a for a of meta).length != 0
        if typeof meta is 'string'
            return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta}\""

        else  # qc
            (k for own k of meta).length isnt 0
            meta_table = '<table class=\'qc_meta_table\'>\n<tr><td></td>'
            for novelty, values of meta when novelty isnt 'all'
                meta_table += "<td>#{novelty}</td>"
            meta_table += '</tr>\n'

            for novelty, val_by_db of meta
                dbs = (db for db, val of val_by_db when db isnt 'average')
                dbs.push 'average'
                break

            short_table = true
            for novelty, val_by_db of meta
                if not check_all_values_equal (val for db, val of val_by_db when db isnt 'average')
                    short_table = false

            if short_table  # Values are the same for each database
                meta_table += '<tr><td></td>'
                for novelty, val_by_db of meta when novelty isnt 'all'
                    pretty_str = toPrettyString(val_by_db[dbs[0]], rec.metric.unit)
                    meta_table += "<td>#{pretty_str}</td>"
                meta_table += '</tr>\n'
            else
                for db in dbs
                    meta_table += "<tr><td>#{db}</td>"
                    for novelty, val_by_db of meta when novelty isnt 'all'
                        meta_table += "<td>#{toPrettyString(val_by_db[db], rec.metric.unit)}</td>"
                    meta_table += '</tr>\n'

            meta_table += '</table>\n'

            return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta_table}\""
    else
        return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\""


get_metric_name_html = (metric, use_full_name=false) ->
    if metric.short_name and not use_full_name
        metricName = metric.short_name
        description = metric.description or metric.name
        return "<a class=\"metric_name\" rel=\"tooltip\" title=\"#{description}\">#{metricName}</a>"
    else
        return metric.name


calc_record_cell_contents = (rec, font) ->
    value = rec.value
    num_html = ''

    if not value? or value == '' or not rec.metric
        rec.cell_contents = '-'

    else
        if typeof value == 'number'
            rec.num = value
            rec.cell_contents = toPrettyString value, rec.metric.unit
            num_html = toPrettyString value

        else if /^-?.?[0-9]/.test value
            result = /([0-9\.]+)(.*)/.exec value
            rec.num = parseFloat result[1]
            rec.cell_contents = toPrettyString(rec.num, rec.metric.unit) + result[2]
            num_html = toPrettyString(rec.num)
        else
            rec.cell_contents = value

    # Max frac width of column
    rec.frac_width = $.fn.intPartTextWidth num_html, font


mean = (a, b) -> (a + b) / 2


calc_cell_contents = (report, section, font) ->
    max_frac_widths_by_metric = {}

    # First round: calculatings max/min integral/fractional widths (for decimal alingment) and max/min values (for heatmaps)
    if not report.type? or report.type == 'FullReport' or report.type == 'SquareSampleReport'
        for sampleReport in (if report.hasOwnProperty('sample_reports') then report.sample_reports else [report])
            for rec in sampleReport.records
                calc_record_cell_contents rec, font

            for rec in sampleReport.records when rec.metric and rec.metric.name of section.metrics_by_name
                if not (rec.metric.name of max_frac_widths_by_metric)
                    max_frac_widths_by_metric[rec.metric.name] = rec.frac_width
                else if rec.frac_width > max_frac_widths_by_metric[rec.metric.name]
                    max_frac_widths_by_metric[rec.metric.name] = rec.frac_width

                rec.metric.values = [] if not rec.metric.values?
                rec.metric.values.push rec.num

    else if report.type == 'SampleReport'
        for rec in report.records when rec.metric.name of section.metrics_by_name
            if rec.num?
                rec.metric.values = [] if not rec.metric.values?
                rec.metric.values.push rec.num

#    else if report.type == 'SquareSampleReport'
#        sampleReport = report
#
#        calc_records_cell_contents sampleReport.records, font
#        for val in  rec in sampleReport.records when rec.metric.name of section.metrics_by_name
#            if not (rec.metric.name of max_frac_widths_by_metric)
#                max_frac_widths_by_metric[rec.metric.name] = rec.frac_width
#            else if rec.frac_width > max_frac_widths_by_metric[rec.metric.name]
#                max_frac_widths_by_metric[rec.metric.name] = rec.frac_width
#
#            if rec.num?
#                rec.metric.values = [] if not rec.metric.values?
#                rec.metric.values.push rec.num


    for metric in section.metrics when metric.values?
        vals = metric.values.slice().sort((a, b) -> if a? and b? then a - b else if a? then a else b)
        l = vals.length

        metric.min = vals[0]
        metric.max = vals[l - 1]
        metric.all_values_equal = metric.min == metric.max
        metric.med = if l % 2 != 0 then vals[(l - 1) / 2] else mean(vals[l / 2], vals[(l / 2) - 1])
        q1 = vals[Math.floor((l - 1) / 4)]
        q3 = vals[Math.floor((l - 1) * 3 / 4)]

        d = q3 - q1
        metric.low_outer_fence = q1 - 3   * d
        metric.low_inner_fence = q1 - 1.5 * d
        metric.top_inner_fence = q3 + 1.5 * d
        metric.top_outer_fence = q3 + 3   * d

    # Second round: setting shift and color properties based on max/min widths and vals
    for sampleReport in (if report.hasOwnProperty('sample_reports') then report.sample_reports else [report])
        for rec in sampleReport.records when rec.metric and rec.metric.name of section.metrics_by_name
            # Padding based on frac width
            if rec.frac_width?
                rec.right_shift = max_frac_widths_by_metric[rec.metric.name] - rec.frac_width

            metric = rec.metric

            # Color heatmap
            if rec.num?
                [top_hue, inner_top_brt, outer_top_brt] = [BLUE_HUE, BLUE_INNER_BRT, BLUE_OUTER_BRT]
                [low_hue, inner_low_brt, outer_low_brt] = [RED_HUE, RED_INNER_BRT, RED_OUTER_BRT]

                if metric.quality == 'Less is better'  # then swap colors
                    [top_hue, low_hue] = [low_hue, top_hue]
                    [inner_top_brt, inner_low_brt] = [inner_low_brt, inner_top_brt]
                    [outer_top_brt, outer_low_brt] = [outer_low_brt, outer_top_brt]

                if not metric.all_values_equal
                    rec.text_color = 'black'

                    # Low outliers
                    if rec.num < rec.metric.low_outer_fence
                        rec.color = get_color low_hue, outer_low_brt
                        rec.text_color = 'white'

                    else if rec.num < rec.metric.low_inner_fence
                        rec.color = get_color low_hue, inner_low_brt

                    # Normal values
                    else if rec.num < metric.med
                        k = (MEDIAN_BRT - MIN_NORMAL_BRT) / (metric.med - rec.metric.low_inner_fence)
                        brt = Math.round MEDIAN_BRT - (metric.med - rec.num) * k
                        rec.color = get_color low_hue, brt

                    # High outliers
                    else if rec.num > rec.metric.top_inner_fence
                        rec.color = get_color top_hue, inner_top_brt

                    else if rec.num > rec.metric.top_outer_fence
                        rec.color = get_color top_hue, outer_top_brt
                        rec.text_color = 'white'

                    else if rec.num > metric.med
                        k = (MEDIAN_BRT - MIN_NORMAL_BRT) / (rec.metric.top_inner_fence - metric.med)
                        brt = Math.round MEDIAN_BRT - (rec.num - metric.med) * k
                        rec.color = get_color top_hue, brt
    return report


median = (x) ->
    return null if (x.length == 0)
    sorted = x.slice().sort((a, b) -> a - b)
    if sorted.length % 2 == 1
        sorted[(sorted.length - 1) / 2]
    else
        (sorted[(sorted.length / 2) - 1] + sorted[(sorted.length / 2)]) / 2


reporting.buildTotalReport = (report, section, columnOrder) ->
    if section.title?
        $('#report').append "<h3 class='table_name'>#{section.title}</h3>"

    calc_cell_contents report, section, $('#report').css 'font'

    table = "<table cellspacing=\"0\"
                    class=\"report_table tableSorter #{if DRAGGABLE_COLUMNS then 'draggable' else ''} fix-align-char\"
                    id=\"report_table_#{section.name}\">"
    table += "\n<thead><tr class=\"top_row_tr\">"
    table += "<th class=\"top_left_td left_column_td\" data-sortBy='numeric'>
                    <span>Sample</span>
              </th>"

    for colNum in [0...section.metrics.length]
        pos = columnOrder[colNum]
        metric = section.metrics[pos]
        if not metric.values?
            continue
        sort_by = if metric.all_values_equal then 'nosort' else 'numeric'
        direction = if metric.quality == 'Less is better' then 'ascending' else 'descending'
        table += "<th class='second_through_last_col_headers_td' data-sortBy=#{sort_by} data-direction=#{direction}position='#{pos}'>
             <span class=\'metricName #{if DRAGGABLE_COLUMNS then 'drag_handle' else ''}\'>#{get_metric_name_html(metric)}</span>
        </th>"
        #{if DRAGGABLE_COLUMNS then '<span class=\'drag_handle\'><span class=\'drag_image\'></span></span>' else ''}
    table += '</tr></thead><tbody>'

    i = 0
    sample_reports_length = if report.hasOwnProperty('sample_reports') then report.sample_reports.length else 1
    for sampleReport in (if report.hasOwnProperty('sample_reports') then report.sample_reports else [report])
        line_caption = sampleReport.display_name  # sample name
        max_sample_name_len = 50
        if line_caption.length > max_sample_name_len
            line_caption = "<span title=\"#{line_caption}\">#{line_caption.substring(0, max_sample_name_len)}...</span>"

        second_row_tr = if i == 0 then "second_row_tr" else ""
        table += "\n<tr class=\"#{second_row_tr}\">
            <td class=\"left_column_td td\" data-sortAs=#{sample_reports_length - i}>"
        if sample_reports_length == 1
            table += "<span class=\"sample_name\">#{line_caption}</span>"
        else
            if sampleReport.html_fpath?
                if typeof sampleReport.html_fpath is 'string'
                    table += "<a class=\"sample_name\" href=\"#{sampleReport.html_fpath}\">#{line_caption}</a>"
                else  # several links for one sample are possible multi-reports (e.g. TargQC)
                    if (k for own k of sampleReport.html_fpath).length == 0
                        table += "<span class=\"sample_name\">#{line_caption}</span>"
                    else
                        links = ""
                        for report_name, html_fpath of sampleReport.html_fpath
                            if not html_fpath
                              continue
                            if links.length != 0
                              links += ", "
                            links += "<a href=\"#{html_fpath}\">#{report_name}</a>"
                        table += "<span class=\"sample_name\">#{line_caption} (#{links})</span>"
                # table += "<a class=\"sample_name\" href=\"#{sampleReport.html_fpath}\">#{line_caption}</a>"
            else
                table += "<span class=\"sample_name\">#{line_caption}</span>"

        table += "</td>"

        for colNum in [0...section.metrics.length]
            pos = columnOrder[colNum]
            metric = section.metrics[pos]
            if not metric.values?
                continue
            rec = null
            for r in sampleReport.records
                if !r.metric
                    continue
                if r.metric.name == metric.name
                    rec = r
                    break
            if not rec?
                table += "<td></td>"
                continue

            table += "<td metric=\"#{metric.name}\"
                          style=\"background-color: #{rec.color}; color: #{rec.text_color}\"
                          quality=\"#{metric.quality}\"
                          class='td "
            if rec.num?
                table += " number' number=\"#{rec.value}\" data-sortAs=#{rec.value}>"
            else
                table += "'>"

            if rec.right_shift?
                padding = "margin-left: #{rec.right_shift}px; margin-right: -#{rec.right_shift}px;"
            else
                padding = ""

            if rec.html_fpath?
                if typeof rec.html_fpath is 'string'
                    table += "<a href=\"#{rec.html_fpath}\">#{rec.cell_contents}
                            </a>
                          </td>"
                else  # varQC -- several variant callers for one sample are possible
                    if (k for own k of rec.html_fpath).length == 0
                        rec.value = null
                        calc_record_cell_contents rec, $('#report').css 'font'
                        table += "#{rec.cell_contents}</td>"
                    else
                        caller_links = ""
                        for caller, html_fpath of rec.html_fpath
                            if caller_links.length != 0
                                caller_links += ", "
                                caller_links += "<a href=\"#{html_fpath}\">#{caller}</a>"
                        table += "#{rec.cell_contents} (#{caller_links})</td>"
            else
                table += "<a style=\"#{padding}\"
                              #{get_meta_tag_contents(rec)}>#{rec.cell_contents}
                            </a>
                          </td>"
        table += "</tr>"
        i += 1
    table += "\n</tbody></table>\n"
#    table += "<div style=\"height: 30px; display: block;\"></div>"

    $('#report').append table


reporting.buildCommonRecords = (common_records) ->
    if common_records.length == 0
        return

    for rec in common_records
        calc_record_cell_contents rec, $('#report').css 'font'

    table = "<table cellspacing=\"0\" class=\"common_table\" id=\"common_table\">"
    for rec in common_records
        table += "\n<tr><td>"
        if rec.html_fpath?
            table += "<a href=\"#{rec.html_fpath}\">
                         #{rec.cell_contents}</a>"
        else
            table += "<span class='metric_name'>#{get_metric_name_html(rec.metric, use_full_name=true)}:</span>
                          #{rec.cell_contents}"
        table += "</td></tr>"
    table += "\n</table>\n"
    table += "<div class='space_30px'></div>"

    $('#report').append table


#set_legend = ->
#    legend = '<span>'
#    step = 6
#    for hue in [RED_HUE..GREEN_HUE] by step
#        legend += "<span style=\"background-color: #{get_color hue}\">"
#
#        switch hue
#            when RED_HUE              then legend += 'w'
#            when RED_HUE   +     step then legend += 'o'
#            when RED_HUE   + 2 * step then legend += 'r'
#            when RED_HUE   + 3 * step then legend += 's'
#            when RED_HUE   + 4 * step then legend += 't'
#            when GREEN_HUE - 3 * step then legend += 'b'
#            when GREEN_HUE - 2 * step then legend += 'e'
#            when GREEN_HUE -     step then legend += 's'
#            when GREEN_HUE            then legend += 't'
#            else                           legend += '.'
#        legend += "</span>"
#    legend += "</span>"
#    $('#report_legend').append legend


$.fn._splitDot_partTextWidth = (html, font, part_type) ->  # part_type = 'int'|'frac'
    parts = html.split '.'

    if part_type == 'frac'
        if parts.length < 2
            return 0
        else
            frac_part = '.' + parts[1]

    else if part_type == 'int'
        frac_part = parts[0]

    if (!$.fn.fracPartTextWidth.fakeEl)
        $.fn.fracPartTextWidth.fakeEl = $('<span>').hide().appendTo document.body

    $.fn.fracPartTextWidth.fakeEl.html frac_part
    $.fn.fracPartTextWidth.fakeEl.css 'font', font
    return $.fn.fracPartTextWidth.fakeEl.width()


$.fn.fracPartTextWidth = (html, font) ->
    $.fn._splitDot_partTextWidth html, font, 'frac'


$.fn.intPartTextWidth = (html, font) ->
    $.fn._splitDot_partTextWidth html, font, 'int'


$.fn.textWidth = (text, font) ->
    if (!$.fn.textWidth.fakeEl)
        $.fn.textWidth.fakeEl = $('<span>').hide().appendTo document.body

    $.fn.textWidth.fakeEl.html text
    $.fn.textWidth.fakeEl.css 'font', font
    return $.fn.textWidth.fakeEl.width()


String.prototype.trunc = (n) ->
    this.substr(0, n - 1) + (this.length > n ? '&hellip;': '')


#postprocess_cells = ->
#    processes_metrics = []
#
#    $(".report_table td[number]").each ->
#        metricName = $(this).attr 'metricName'
#
#        if !(metricName in processes_metrics)
#            processes_metrics.push metricName
#            console.log metricName
#
#            quality = $(this).attr 'quality'
#            all_cells = $('.report_table').find "td[metricName=\"#{metricName}\"][number]"
#            all_numbers = ($(cell).attr 'number' for cell in all_cells)
#
#            set_heatmap all_cells, all_numbers, quality
#
#            set_offset all_cells, all_numbers, metricName