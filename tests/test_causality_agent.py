import json
from kqml import KQMLList, KQMLString, KQMLPerformative
from indra.statements import stmts_from_json
from causality_agent.causality_module import _resource_dir
from causality_agent import causality_agent
from causality_agent.causality_module import CausalityModule
from bioagents.tests.integration import _IntegrationTest
from bioagents.tests.util import ekb_kstring_from_text, ekb_from_text, get_request, agent_clj_from_text
import time

ca = causality_agent.CausalityAgent(_resource_dir)

def _get_names_from_list_cljson(list_cljson):
    return list(map(lambda o: o.gets('NAME'), list_cljson))

def _read_from_kqml_list(kl, prop_path):
    """read a property of a kqml list from given path array"""
    res = kl
    for p in prop_path:
        res = res.get(p)
    return res

def _reads_from_kqml_list(kl, prop_path):
    """read a property of a kqml list from given path array as string"""
    val = _read_from_kqml_list(kl, prop_path)
    if not val:
        return None
    return val.string_value()

class TestCausalPath(_IntegrationTest):
    def __init__(self, *args):
        super(TestCausalPath, self).__init__(CausalityModule)

    def create_message(self):
        source = agent_clj_from_text('MAPK1')
        target = agent_clj_from_text('JUND')
        content = KQMLList('FIND-CAUSAL-PATH')
        content.set('source', source)
        content.set('target', target)
        content.sets('direction', 'both')

        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS'
        paths = output.get('paths')
        assert len(paths) == 1
        path = paths[0]
        assert _reads_from_kqml_list(path, ['enz', 'name']) == 'MAPK1'
        assert _reads_from_kqml_list(path, ['sub', 'name']) == 'JUND'


    def create_message_2(self):
        source = agent_clj_from_text('MAPK1')
        target = agent_clj_from_text('CREB1')
        content = KQMLList('FIND-CAUSAL-PATH')
        content.set('source', source)
        content.set('target', target)
        content.sets('direction', 'both')
        msg = get_request(content)
        return msg, content

    def check_response_to_message_2(self, output):
        assert output.head() == 'SUCCESS'
        paths = output.get('paths')
        path = paths[0]
        # assert len(paths) == 1
        assert _reads_from_kqml_list(path, ['enz', 'name']) == 'MAPK1'
        assert _reads_from_kqml_list(path, ['sub', 'name']) == 'CREB1'
        assert _reads_from_kqml_list(path, ['residue']) == 'S'
        assert _reads_from_kqml_list(path, ['position']) == '133'


    def create_message_failure(self):
        source = agent_clj_from_text('MAPK1')
        target = agent_clj_from_text('RAS')
        content = KQMLList('FIND-CAUSAL-PATH')
        content.set('source', source)
        content.set('target', target)
        content.sets('direction', 'both')
        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure(self, output):
        assert output.head() == 'FAILURE'
        reason = output.gets('reason')
        assert reason == 'NO_PATH_FOUND'

    # def create_message_failure_2(self):
    #     source = ekb_kstring_from_text('ABC')
    #     target = ekb_kstring_from_text('TP53')
    #     content = KQMLList('FIND-CAUSAL-PATH')
    #     content.set('source', source)
    #     content.set('target', target)
    #     content.sets('direction', 'both')
    #     msg = get_request(content)
    #     return msg, content
    #
    # def check_response_to_message_failure_2(self, output):
    #     assert output.head() == 'FAILURE'
    #     reason = output.gets('reason')
    #     assert reason == 'NO_PATH_FOUND'

    def create_message_failure_3(self):
        source = agent_clj_from_text('RAS')
        target = agent_clj_from_text('MAPK1')
        content = KQMLList('FIND-CAUSAL-PATH')
        content.set('source', source)
        content.set('target', target)
        content.sets('direction', 'strict')
        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure_3(self, output):
        assert output.head() == 'FAILURE'
        reason = output.gets('reason')
        assert reason == 'NO_PATH_FOUND'


class TestCausalityTarget(_IntegrationTest):
    def __init__(self, *args):
        super(TestCausalityTarget, self).__init__(CausalityModule)

    def create_message(self):
        source = agent_clj_from_text('MAPK1')
        content = KQMLList('FIND-CAUSALITY-TARGET')
        content.set('source', source)
        content.sets('type', 'phosphorylation')
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        paths = output.get('paths')
        assert len(paths) == 904
        path = paths[0]
        assert _reads_from_kqml_list(path, ['sub', 'name']) == 'RPTOR'
        assert _reads_from_kqml_list(path, ['residue']) == 'S'
        assert _reads_from_kqml_list(path, ['position']) == '863'

    def create_message_failure(self):
        source = agent_clj_from_text('MAPK1')
        content = KQMLList('FIND-CAUSALITY-TARGET')
        content.set('source', source)
        content.sets('type', 'activation')
        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure(self, output):
        assert output.head() == 'FAILURE', output
        reason = output.gets('reason')
        assert reason == "MISSING_MECHANISM"

    # def create_message_failure_2(self):
    #     target = ekb_kstring_from_text('ABC')
    #     content = KQMLList('FIND-CAUSALITY-TARGET')
    #     content.set('target', target)
    #     content.sets('type', 'phosphorylates')
    #     msg = get_request(content)
    #     return msg, content
    #
    # def check_response_to_message_failure_2(self, output):
    #     assert output.head() == 'FAILURE', output
    #     reason = output.gets('reason')
    #     assert reason == "MISSING_MECHANISM"


class TestCausalitySource(_IntegrationTest):
    def __init__(self, *args):
        super(TestCausalitySource, self).__init__(CausalityModule)

    def create_message(self):
        target = agent_clj_from_text('BRAF')
        content = KQMLList('FIND-CAUSALITY-SOURCE')
        content.set('target', target)
        content.sets('type', 'phosphorylation')
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        paths = output.get('paths')
        assert len(paths) == 80
        path = paths[0]
        assert _reads_from_kqml_list(path, ['enz', 'name']) == 'NRAS'
        assert _reads_from_kqml_list(path, ['residue']) == 'S'
        assert _reads_from_kqml_list(path, ['position']) == '601'

    def create_message_failure(self):
        target = agent_clj_from_text('BRAF')
        content = KQMLList('FIND-CAUSALITY-SOURCE')
        content.set('target', target)
        content.sets('type', 'activates')
        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure(self, output):
        assert output.head() == 'FAILURE', output
        reason = output.gets('reason')
        assert reason == "MISSING_MECHANISM"

    # def create_message_failure_2(self):
    #     target = ekb_kstring_from_text('ABC')
    #     content = KQMLList('FIND-CAUSALITY-SOURCE')
    #     content.set('target', target)
    #     content.sets('type', 'phosphorylates')
    #     msg = get_request(content)
    #     return msg, content
    #
    # def check_response_to_message_failure_2(self, output):
    #     assert output.head() == 'FAILURE', output
    #     reason = output.gets('reason')
    #     assert reason == "MISSING_MECHANISM"


class TestNextCorrelation(_IntegrationTest):
    def __init__(self, *args):
        super(TestNextCorrelation, self).__init__(CausalityModule)

    def create_message_01_explainable(self):
        source = agent_clj_from_text('AKT1')
        content = KQMLList('DATASET-CORRELATED-ENTITY')
        content.set('source', source)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_01_explainable(self, output):
        assert output.head() == 'SUCCESS', output
        target = output.gets('target')
        correlation = output.gets('correlation')
        explainable = output.gets('explainable')


        assert target == 'BRAF'
        assert correlation == str(0.7610843243760473)
        assert explainable == 'explainable'


    def create_message_02_explainable2(self):
        time.sleep(2)
        source = agent_clj_from_text('AKT1')
        content = KQMLList('DATASET-CORRELATED-ENTITY')
        content.set('source', source)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_02_explainable2(self, output):
        assert output.head() == 'SUCCESS', output
        target = output.gets('target')
        correlation = output.gets('correlation')
        explainable = output.gets('explainable')

        assert target == 'PTPN1'
        assert correlation.startswith('0.581061418186')
        assert explainable == 'explainable'


    def create_message_03_unexplainable(self):
        time.sleep(2)
        source = agent_clj_from_text('AKT1')
        content = KQMLList('DATASET-CORRELATED-ENTITY')
        content.set('source', source)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_03_unexplainable(self, output):
        assert output.head() == 'SUCCESS', output
        target = output.gets('target')
        correlation = output.gets('correlation')
        explainable = output.gets('explainable')
        assert target == 'AGPS'
        assert correlation.startswith('0.94999636806')
        assert explainable == 'unexplainable'
        time.sleep(1)

    def create_message_04_reset(self):
        time.sleep(2)
        content = KQMLList('RESET-CAUSALITY-INDICES')
        msg = get_request(content)
        return msg, content

    def check_response_to_message_04_reset(self, output):
        assert output.head() == 'SUCCESS', output


    def create_message_05_explainable_again(self):
        time.sleep(2)
        source = agent_clj_from_text('AKT1')
        content = KQMLList('DATASET-CORRELATED-ENTITY')
        content.set('source', source)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_05_explainable_again(self, output):
        assert output.head() == 'SUCCESS', output
        target = output.gets('target')
        correlation = output.gets('correlation')
        explainable = output.gets('explainable')
        assert target == 'BRAF'
        assert correlation.startswith('0.7610843243760473')
        assert explainable == 'explainable'

    # def create_message_failure(self):
    #     time.sleep(2)
    #     source = ekb_kstring_from_text('ABC')
    #     content = KQMLList('DATASET-CORRELATED-ENTITY')
    #     content.set('source', source)
    #     msg = get_request(content)
    #     return msg, content
    #
    # def check_response_to_message_failure(self, output):
    #     assert output.head() == 'FAILURE', output
    #     reason = output.gets('reason')
    #     assert reason == "NO_PATH_FOUND"




class TestCommonUpstreams(_IntegrationTest):
    def __init__(self, *args):
        super(TestCommonUpstreams, self).__init__(CausalityModule)

    def create_message(self):
        content = KQMLList('FIND-COMMON-UPSTREAMS')
        genes = KQMLList([agent_clj_from_text('AKT1'), agent_clj_from_text('BRAF'), agent_clj_from_text('MAPK1')])
        content.set('genes', genes)
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        upstreams = output.get('upstreams')
        assert 'EGF' in _get_names_from_list_cljson(upstreams)

    def create_message_failure(self):
        content = KQMLList('FIND-COMMON-UPSTREAMS')
        genes = KQMLList([agent_clj_from_text('UGT2B10'), agent_clj_from_text('PTEN')])
        content.set('genes', genes)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure(self, output):
        assert output.head() == 'FAILURE', output
        reason = output.gets('reason')
        assert reason == "NO_UPSTREAM_FOUND"


class TestMutex(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutex, self).__init__(CausalityModule)

    def create_message(self):
        content = KQMLList('FIND-MUTEX')
        gene = agent_clj_from_text('TP53')
        disease = agent_clj_from_text('breast cancer')
        content.set('gene', gene)
        content.set('disease', disease)
        msg = get_request(content)
        return msg, content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        mutex = output.get('mutex')


        test_res = KQMLList.from_string('((:score 0.0 :group (TP53 CDH1)) '
                                        '(:score 0.0 :group (CDH1 TP53)) '
                                        '(:score 0.0 :group (GATA3 TP53 CDH1)) '
                                        '(:score 0.0 :group (CTCF TP53 CDH1 GATA3)))')

        # TODO: do this without converting into string
        assert str(test_res) == str(mutex)

    def create_message_failure(self):
        content = KQMLList('FIND-MUTEX')
        gene = agent_clj_from_text('BRAF')
        content.set('gene', gene)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure(self, output):
        assert output.head() == 'FAILURE', output
        reason = output.gets('reason')
        assert reason == "MISSING_MECHANISM"


class TestMutSig(_IntegrationTest):
    def __init__(self, *args):
        super(TestMutSig, self).__init__(CausalityModule)

    def create_message_OV(self):
        content = KQMLList('FIND-MUTATION-SIGNIFICANCE')
        gene = agent_clj_from_text('TP53')
        disease = agent_clj_from_text('Ovarian serous cystadenocarcinoma')
        content.set('gene', gene)
        content.set('disease', disease)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_OV(self, output):
        assert output.head() == 'SUCCESS', output
        mut_sig = output.gets('mutsig')
        assert mut_sig == "highly significant"

    def create_message_PAAD(self):
        content = KQMLList('FIND-MUTATION-SIGNIFICANCE')
        gene = agent_clj_from_text('ACTN4')
        disease = agent_clj_from_text('pancreatic adenocarcinoma')
        content.set('gene', gene)
        content.set('disease', disease)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_PAAD(self, output):
        assert output.head() == 'SUCCESS', output
        mut_sig = output.gets('mutsig')
        assert mut_sig == "not significant"

    def create_message_failure(self):
        content = KQMLList('FIND-MUTATION-SIGNIFICANCE')
        gene = agent_clj_from_text('ACTN4')
        disease = agent_clj_from_text('abc cancer')
        content.set('gene', gene)
        content.set('disease', disease)
        msg = get_request(content)
        return msg, content

    def check_response_to_message_failure(self, output):
        assert output.head() == 'FAILURE', output
        reason = output.gets('reason')
        assert reason == "INVALID_DISEASE"

    # def create_message_failure_2(self):
    #     content = KQMLList('FIND-MUTATION-SIGNIFICANCE')
    #     gene = ekb_kstring_from_text('ABC')
    #     disease = ekb_from_text('breast cancer')
    #     content.set('gene', gene)
    #     content.set('disease', disease)
    #     msg = get_request(content)
    #     return msg, content
    #
    # def check_response_to_message_failure_2(self, output):
    #     assert output.head() == 'FAILURE', output
    #     reason = output.gets('reason')
    #     assert reason == "MISSING_MECHANISM"

class TestCellularLocation(_IntegrationTest):
    def __init__(self, *args):
        super(TestCellularLocation, self).__init__(CausalityModule)

    def create_message_AKT1(self):
        content = KQMLList('FIND-CELLULAR-LOCATION')
        genes = agent_clj_from_text('AKT1')
        content.set('genes', genes)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_AKT1(self, output):
        assert output.head() == 'SUCCESS', output
        components = output.get('components')
        assert 'mitochondrion' in _get_names_from_list_cljson(components)

    def create_message_2(self):
        content = KQMLList('FIND-CELLULAR-LOCATION-FROM-NAMES')
        content.set('genes', ['AKT1', 'MAPK1'])

        msg = get_request(content)
        return msg, content

    def check_response_to_message_2(self, output):
        assert output.head() == 'SUCCESS', output
        components = output.get('components')
        assert 'mitochondrion' in _get_names_from_list_cljson(components)

    def create_message_3(self):
        content = KQMLList('FIND-CELLULAR-LOCATION')
        genes = KQMLList([agent_clj_from_text('AKT1'), agent_clj_from_text('MAPK1')])
        content.set('genes', genes)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_3(self, output):
        assert output.head() == 'SUCCESS', output
        components = output.get('components')
        assert 'mitochondrion' in _get_names_from_list_cljson(components)

    def create_message_4(self):
        content = KQMLList('FIND-CELLULAR-LOCATION-FROM-NAMES')
        content.set('genes', ['MAPK1',  'RAS'])

        msg = get_request(content)
        return msg, content

    def check_response_to_message_4(self, output):
        assert output.head() == 'FAILURE', output
        reason = output.gets('reason')
        assert reason == "NO_COMMON_CELLULAR_LOCATION_FOUND"

class TestGeneSummary(_IntegrationTest):

    def __init__(self, *args):
        super(TestGeneSummary, self).__init__(CausalityModule)

    def create_message_AKT1(self):
        content = KQMLList('FIND-GENE-SUMMARY')
        genes = agent_clj_from_text('AKT1')
        content.set('gene', genes)

        msg = get_request(content)
        return msg, content

    def check_response_to_message_AKT1(self, output):
        #TODO: failing because of 503 error
        assert output.head() == 'SUCCESS', output
        components = output.get('geneSummary')
        assert 'AKT1' in components.data
