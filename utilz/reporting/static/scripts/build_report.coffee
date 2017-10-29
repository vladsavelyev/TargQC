showPlotWithInfo = (info) ->
    newSeries = []
    newColors = []

    $('#legend-placeholder').find 'input:checked'.each ->
        number = $(this).attr 'name'
        if number and info.series && info.series.length > 0
            for i in [i...info.series.length]
                series = info.series[i]
                break if series.number != number

            if i <= info.series.length
                newSeries.push(series)
                newColors.push(series.color)
            else
                if window.console? then console.log('no series with number ' + number)

        if newSeries.length == 0
            newSeries.push
                data: []

            newColors.push '#FFF'

        info.showWithData(newSeries, newColors)


recoverOrderFromCookies = (report_name) ->
    return null unless navigator.cookieEnabled

    orderString = readCookie report_name + '_order'
    return null unless orderString

    columnOrder = []
    fail = false

    for val in orderString.split(' ')
        val = parseInt val
        if isNaN val
            fail = true
        else
            columnOrder.push val

    return null if fail

    return columnOrder


readJson = (what) ->
    result
    try
        result = JSON.parse $('#' + what + '-json').html()
    catch e
        result = null

    return result


totalReportData =
    date: null
    report: null

report =
    name: ''
    order: null
    sample_reports: []
    metric_storage:
        general_section:
            name: ''
            metrics: []
        sections: []

section =
    name: ''
    title: ''
    metrics: []
    metrics_by_name: {}

metric =
    name: ''
    short_name: ''
    description: ''
    quality: ''
    common: true
    unit: ''


extend = (object, properties) ->
    for key, val of properties
        object[key] = val
    object


merge = (options, overrides) ->
    extend (extend {}, options), overrides


preprocessReport = (report) ->
    all_metrics_by_name = {}
    for m in report.metric_storage.general_section.metrics
        report.metric_storage.general_section.metrics_by_name[m.name] = m
    extend all_metrics_by_name, report.metric_storage.general_section.metrics_by_name

    for s in report.metric_storage.sections
        for m in s.metrics
            s.metrics_by_name[m.name] = m
        extend all_metrics_by_name, s.metrics_by_name

    if report.hasOwnProperty('sample_reports')
        for sample_report in report.sample_reports
            sample_report.metric_storage = report.metric_storage
            for rec in sample_report.records
                rec.metric = all_metrics_by_name[rec.metric.name]
    else
        for rec in report.records
            rec.metric = all_metrics_by_name[rec.metric.name]

    return report


reporting.buildReport = ->
    totalReportData = readJson 'total-report'
    unless (totalReportData)
        if window.console? then console.log "Error: cannot read #total-report-json"
        return 1

    report = preprocessReport totalReportData.report

    $('#report_date').html '<p>' + totalReportData.date + '</p>'

    common_metrics_by_name = report.metric_storage.general_section.metrics_by_name
    records = if report.hasOwnProperty('sample_reports') then report.sample_reports[0].records else report.records
    general_records = (rec for rec in records when rec and rec.metric and rec.metric.name of common_metrics_by_name)
    reporting.buildCommonRecords general_records

    for section in report.metric_storage.sections
        columnNames = (m.name for m in section.metrics)
        columnOrder = (recoverOrderFromCookies section.name) or report.order or [0...columnNames.length]

        reporting.buildTotalReport report, section, columnOrder
        if report.hasOwnProperty('plots')
            plots_html = ""
            for plot in report.plots
                plots_html += "<img src=\"#{plot}\"/>"
            $('#plot').html plots_html

    return 0