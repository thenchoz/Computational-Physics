% This script tests eig_power.m for correctness

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

test_case_ans = { 2, 3, 10, 2, 3, 10, 10.1, -9.9 };

for i=1:size(test_case_names)
    
    fprintf('Test %d: %s ', i, test_case_names{i});
    
    rng(i);
    
    try
        [~, m_rl, m_cp] = evalc(test_case_names{i});
    catch e
        [~, m_rl] = evalc(test_case_names{i});
        m_cp = m_rl;
    end
    
    [vec, val] = eig_power(m_cp);
    delta = m_cp*vec - val*vec;
    
    assert(abs(val-test_case_ans{i})<1e-10);
    assert(all(abs(delta) < 1e-5));
    
    fprintf('OK\n');
end