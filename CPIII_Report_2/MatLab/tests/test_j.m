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

for i=1:size(test_case_names)
    
    fprintf('Test %d: %s ', i, test_case_names{i});
    
    rng(i);
    
    try
        [~, m_rl, m_cp] = evalc(test_case_names{i});
    catch e
        [~, m_rl] = evalc(test_case_names{i});
        m_cp = m_rl;
    end
    reference = sort(eig(m_rl));
    
    val = sort(eig_j(m_rl));
    
    assert(all(abs(val-reference)<1e-5));
    
    fprintf('OK\n');
end