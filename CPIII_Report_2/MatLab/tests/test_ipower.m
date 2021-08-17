% This script tests eig_ipower.m for correctness

test_case_names = {
    'diag(1:2)';
    'diag(1:3)';
    'diag(1:10)';
    'rmg(1:2)';
    'rmg(1:3)';
    'rmg(1:10)';
    'rmg(-2.9:10.1)';
    'rmg(-9.9:3.1)';
};

test_case_inp = { 1.2, 1.8, 5.6, 1.2, 2.8, -3, 0, 0};
test_case_ans = { 1.0, 2.0, 6.0, 1.0, 3.0, 1.0, 0.1, 0.1 };

for i=1:size(test_case_names)
    
    fprintf('Test %d: %s ', i, test_case_names{i});
    
    rng(i);
    
    try
        [~, m_rl, m_cp] = evalc(test_case_names{i});
    catch e
        [~, m_rl] = evalc(test_case_names{i});
        m_cp = m_rl;
    end
    
    [vec, val] = eig_ipower(m_cp,test_case_inp{i});
    delta = m_cp*vec - val*vec;
    
    assert(abs(val-test_case_ans{i})<1e-10);
    assert(all(abs(delta) < 1e-5));
    
    fprintf('OK\n');
end