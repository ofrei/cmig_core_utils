function Manufacturer = NQLQ_get_manufacturer(str)

if ~isempty(regexpi(str,'ge medical')), Manufacturer = 'ge medical';
  elseif ~isempty(regexpi(str,'siemens')), Manufacturer = 'siemens';
  elseif ~isempty(regexpi(str,'philips')), Manufacturer = 'philips';
  else, Manufacturer = 'unknown';
end

