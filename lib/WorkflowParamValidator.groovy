import java.util.regex.Pattern

class WorkflowParamValidator {
    private static final Pattern TOKEN = Pattern.compile(/[A-Za-z0-9][A-Za-z0-9._+-]*/)
    private static final Pattern INTEGER = Pattern.compile(/[0-9]+/)
    private static final Pattern DECIMAL = Pattern.compile(/[0-9]+(\.[0-9]+)?/)
    private static final Pattern PATH_VALUE = Pattern.compile(/[A-Za-z0-9._+,:=@%\/-]+/)
    private static final Pattern HEADER_VALUE = Pattern.compile(/[^\p{Cntrl}$`"';|&<>]+/)
    private static final Set BOOLEAN_STRINGS = ['TRUE', 'FALSE'] as Set

    static void validate(def params) {
        requirePath(params, 'sdrf')
        requirePath(params, 'resultsRoot')
        requirePath(params, 'referenceFasta')
        requirePath(params, 'transcriptToGene')
        requirePath(params, 'transcriptomeIndex')
        requireToken(params, 'protocol')
        requireToken(params, 'experimentType')

        optionalPath(params, 'manualDownloadFolder')
        optionalPath(params, 'fastqProviderConfig')

        requireToken(params, 'downloadMethod')
        requireInteger(params, 'maxConcurrentDownloads')
        requireInteger(params, 'minMappingRate')
        requireInteger(params, 'minCbFreq')

        requireFields(params.fields, [
            'run',
            'cdna_uri',
            'cell_barcode_uri',
            'cell_barcode_size',
            'umi_barcode_size',
            'end',
            'cell_count'
        ])
        optionalFields(params.fields, ['quality', 'controlled_access', 'techrep'])

        requireNestedProtocol(params, params.protocol.toString())
        validateEmptyDrops(params.emptyDrops)

        if (has(params, 'salmon')) {
            requireNestedInteger(params.salmon?.index, 'salmon.index', 'kmerSize')
        }
    }

    static String safeToken(value, String fieldName) {
        def text = value == null ? '' : value.toString()
        if (!(TOKEN.matcher(text).matches())) {
            throw new IllegalArgumentException("Unsafe SDRF value for ${fieldName}: '${text}'")
        }
        text
    }

    static String safeInteger(value, String fieldName) {
        def text = value == null ? '' : value.toString()
        if (!(INTEGER.matcher(text).matches())) {
            throw new IllegalArgumentException("Unsafe SDRF numeric value for ${fieldName}: '${text}'")
        }
        text
    }

    static String safeUri(value, String fieldName) {
        def text = value == null ? '' : value.toString()
        if (!(text ==~ /[^\p{Cntrl}\s]+/)) {
            throw new IllegalArgumentException("Unsafe SDRF URI value for ${fieldName}: '${text}'")
        }
        text
    }

    static String safeControlledAccess(value) {
        def text = value == null ? 'no' : value.toString().toLowerCase()
        if (!(text in ['yes', 'no'])) {
            throw new IllegalArgumentException("Unsafe SDRF controlled access value: '${value}'")
        }
        text
    }

    static String shellQuote(value) {
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }

    private static void requireNestedProtocol(def params, String protocol) {
        if (!has(params, protocol)) {
            throw new IllegalArgumentException("Missing workflow params.${protocol} protocol settings")
        }
        def protocolParams = params.get(protocol)
        requireNestedToken(protocolParams, "params.${protocol}", 'alevinType')
        requireNestedInteger(protocolParams, "params.${protocol}", 'barcodeLength')
        requireNestedInteger(protocolParams, "params.${protocol}", 'umiLength')
        requireNestedInteger(protocolParams, "params.${protocol}", 'end')
        requireNestedEnum(protocolParams, "params.${protocol}", 'libType', ['ISR', 'ISF', 'IU'] as Set)
        optionalNestedPath(protocolParams, "params.${protocol}", 'whitelist')
    }

    private static void validateEmptyDrops(def emptyDrops) {
        if (emptyDrops == null) {
            throw new IllegalArgumentException('Missing workflow params.emptyDrops settings')
        }
        requireNestedInteger(emptyDrops, 'params.emptyDrops', 'lower')
        requireNestedInteger(emptyDrops, 'params.emptyDrops', 'nIters')
        requireNestedEnum(emptyDrops, 'params.emptyDrops', 'testAmbient', BOOLEAN_STRINGS)
        requireNestedEnum(emptyDrops, 'params.emptyDrops', 'filterEmpty', BOOLEAN_STRINGS)
        requireNestedDecimal(emptyDrops, 'params.emptyDrops', 'filterFdr')
        requireNestedEnum(emptyDrops, 'params.emptyDrops', 'libType', ['ISR', 'ISF', 'IU'] as Set)
    }

    private static void requireFields(def fields, List required) {
        if (fields == null) {
            throw new IllegalArgumentException('Missing workflow params.fields settings')
        }
        required.each { requireNestedHeader(fields, 'params.fields', it) }
    }

    private static void optionalFields(def fields, List optional) {
        optional.each { optionalNestedHeader(fields, 'params.fields', it) }
    }

    private static void requirePath(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", PATH_VALUE)
    }

    private static void optionalPath(def params, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "params.${name}", PATH_VALUE)
        }
    }

    private static void requireToken(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", TOKEN)
    }

    private static void requireInteger(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", INTEGER)
    }

    private static void requireNestedToken(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", TOKEN)
    }

    private static void requireNestedInteger(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", INTEGER)
    }

    private static void requireNestedDecimal(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", DECIMAL)
    }

    private static void requireNestedEnum(def params, String scope, String name, Set allowed) {
        requireValue(params, "${scope}.${name}", name)
        def text = params.get(name).toString()
        if (!(text in allowed)) {
            throw new IllegalArgumentException("${scope}.${name} must be one of ${allowed}; got '${text}'")
        }
    }

    private static void optionalNestedPath(def params, String scope, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "${scope}.${name}", PATH_VALUE)
        }
    }

    private static void requireNestedHeader(def params, String scope, String name) {
        requireValue(params, "${scope}.${name}", name)
        assertPattern(params.get(name), "${scope}.${name}", HEADER_VALUE)
    }

    private static void optionalNestedHeader(def params, String scope, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "${scope}.${name}", HEADER_VALUE)
        }
    }

    private static void requireValue(def params, String label, String key) {
        if (params == null || !has(params, key) || params.get(key) == null || params.get(key).toString() == '') {
            throw new IllegalArgumentException("Missing required workflow parameter ${label}")
        }
    }

    private static void assertPattern(value, String label, Pattern pattern) {
        def text = value == null ? '' : value.toString()
        if (!pattern.matcher(text).matches()) {
            throw new IllegalArgumentException("Invalid workflow parameter ${label}: '${text}'")
        }
    }

    private static boolean has(def params, String key) {
        params != null && params.containsKey(key)
    }
}
