import java.net.InetAddress
import java.net.URI
import java.util.regex.Pattern

class WorkflowParamValidator {
    private static final Pattern TOKEN = Pattern.compile(/[A-Za-z0-9][A-Za-z0-9._+-]*/)
    private static final Pattern INTEGER = Pattern.compile(/[0-9]+/)
    private static final Pattern DECIMAL = Pattern.compile(/[0-9]+(\.[0-9]+)?/)
    private static final Pattern PATH_VALUE = Pattern.compile(/[A-Za-z0-9._+,:=@%\/-]+/)
    private static final Pattern HEADER_VALUE = Pattern.compile(/[^\p{Cntrl}$`"';|&<>]+/)
    private static final Pattern FASTQ_URI_VALUE = Pattern.compile(/[^\p{Cntrl}\s]+/)
    private static final Set DEFAULT_FASTQ_URI_SCHEMES = ['http', 'https', 'ftp', 'sra/http', 'sra/https', 'sra/ftp'] as Set
    private static final Set BOOLEAN_STRINGS = ['TRUE', 'FALSE'] as Set
    private static Set allowedFastqUriSchemes = DEFAULT_FASTQ_URI_SCHEMES
    private static Set allowedFastqUriHosts = [] as Set
    private static boolean denyPrivateFastqNetwork = true

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
        validateFastqDownloadPolicy(params)
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

    static String safeUri(value, String fieldName, String controlledAccess = 'no') {
        def text = value == null ? '' : value.toString()
        if (!(FASTQ_URI_VALUE.matcher(text).matches())) {
            throw new IllegalArgumentException("Unsafe SDRF URI value for ${fieldName}: '${text}'")
        }
        if (controlledAccess != 'yes') {
            validateFastqUriEgress(text, fieldName)
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

    private static void validateFastqDownloadPolicy(def params) {
        def policy = has(params, 'fastqDownload') ? params.fastqDownload : null
        allowedFastqUriSchemes = stringSet(policy != null && has(policy, 'allowedSchemes') ? policy.allowedSchemes : DEFAULT_FASTQ_URI_SCHEMES)
        allowedFastqUriHosts = stringSet(policy != null && has(policy, 'allowedHosts') ? policy.allowedHosts : [])
        denyPrivateFastqNetwork = !(policy != null && has(policy, 'denyPrivateNetwork') && policy.denyPrivateNetwork.toString().toLowerCase() == 'false')

        if (allowedFastqUriSchemes.isEmpty()) {
            throw new IllegalArgumentException('params.fastqDownload.allowedSchemes must contain at least one scheme')
        }
        allowedFastqUriSchemes.each { scheme ->
            if (!(scheme ==~ /(sra\/)?[a-z][a-z0-9+.-]*/)) {
                throw new IllegalArgumentException("Invalid FASTQ download scheme policy entry: '${scheme}'")
            }
        }
        allowedFastqUriHosts.each { host ->
            if (!(host == '*' || host ==~ /(\*\.)?[A-Za-z0-9.-]+/)) {
                throw new IllegalArgumentException("Invalid FASTQ download host policy entry: '${host}'")
            }
        }
    }

    private static void validateFastqUriEgress(String uriValue, String fieldName) {
        def parsed = parseFastqUri(uriValue, fieldName)
        if (!(parsed.scheme in allowedFastqUriSchemes)) {
            throw new IllegalArgumentException("FASTQ URI for ${fieldName} uses disallowed scheme '${parsed.scheme}'")
        }
        if (!hostAllowed(parsed.host)) {
            throw new IllegalArgumentException("FASTQ URI for ${fieldName} targets host '${parsed.host}', which is not allowed by params.fastqDownload.allowedHosts")
        }
        if (denyPrivateFastqNetwork && privateOrLocalHost(parsed.host)) {
            throw new IllegalArgumentException("FASTQ URI for ${fieldName} targets private or local host '${parsed.host}'")
        }
    }

    private static Map parseFastqUri(String uriValue, String fieldName) {
        def candidate = uriValue
        def wrapper = ''
        if (candidate.toLowerCase().startsWith('sra/')) {
            wrapper = 'sra/'
            candidate = candidate.substring(4)
        }
        URI parsed
        try {
            parsed = new URI(candidate)
        } catch (Exception e) {
            throw new IllegalArgumentException("FASTQ URI for ${fieldName} is not parseable: '${uriValue}'")
        }
        def scheme = parsed.scheme == null ? '' : parsed.scheme.toLowerCase()
        def host = parsed.host == null ? '' : parsed.host.toLowerCase().replaceFirst(/\.$/, '')
        if (scheme == '' || host == '') {
            throw new IllegalArgumentException("FASTQ URI for ${fieldName} must include a network scheme and host: '${uriValue}'")
        }
        [scheme: "${wrapper}${scheme}", host: host]
    }

    private static boolean hostAllowed(String host) {
        allowedFastqUriHosts.isEmpty() || allowedFastqUriHosts.any { allowed ->
            allowed == '*' ||
                host == allowed ||
                (allowed.startsWith('*.') && host.endsWith(allowed.substring(1)) && host != allowed.substring(2))
        }
    }

    private static boolean privateOrLocalHost(String host) {
        def h = host.toLowerCase()
        if (h in ['localhost', 'metadata.google.internal']) {
            return true
        }
        if (h.contains(':')) {
            return isPrivateIpv6(h)
        }
        if (!h.contains('.') || h.endsWith('.localhost') || h.endsWith('.local') || h.endsWith('.internal') || h.endsWith('.lan') || h.endsWith('.home') || h.endsWith('.corp') || h.endsWith('.localdomain')) {
            return true
        }
        isPrivateIpv4(h) || isPrivateIpv6(h)
    }

    private static boolean isPrivateIpv4(String host) {
        if (!(host ==~ /[0-9]{1,3}(\.[0-9]{1,3}){3}/)) {
            return false
        }
        def parts = host.split(/\./).collect { it.toInteger() }
        if (parts.any { it < 0 || it > 255 }) {
            return true
        }
        def a = parts[0]
        def b = parts[1]
        def c = parts[2]
        a == 0 ||
            a == 10 ||
            a == 127 ||
            (a == 100 && b >= 64 && b <= 127) ||
            (a == 169 && b == 254) ||
            (a == 172 && b >= 16 && b <= 31) ||
            (a == 192 && b == 168) ||
            (a == 192 && b == 0 && c == 0) ||
            (a == 192 && b == 0 && c == 2) ||
            (a == 198 && b >= 18 && b <= 19) ||
            (a == 198 && b == 51 && c == 100) ||
            (a == 203 && b == 0 && c == 113) ||
            a >= 224
    }

    private static boolean isPrivateIpv6(String host) {
        if (!host.contains(':')) {
            return false
        }
        try {
            def address = InetAddress.getByName(host)
            def lower = host.toLowerCase()
            address.anyLocalAddress ||
                address.loopbackAddress ||
                address.linkLocalAddress ||
                address.siteLocalAddress ||
                address.multicastAddress ||
                lower.startsWith('fc') ||
                lower.startsWith('fd')
        } catch (Exception e) {
            true
        }
    }

    private static Set stringSet(value) {
        def values = value instanceof Collection ? value : value.toString().split(',')
        def result = [] as Set
        values.each { item ->
            def text = item == null ? '' : item.toString().trim().toLowerCase()
            if (text != '') {
                result << text
            }
        }
        result
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
