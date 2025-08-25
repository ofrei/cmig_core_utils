function res = isinteractive()

res = ~isdeployed & isempty(java.lang.System.getProperty( 'java.awt.headless' ));

